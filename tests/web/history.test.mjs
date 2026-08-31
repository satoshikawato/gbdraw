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
  'runtime-test-hooks.js',
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
  let setting = 0;
  let checkpointBuilds = 0;
  const history = createHistoryManager({
    buildIntent: async () => ({ setting }),
    applyIntent: async (intent) => {
      setting = intent.setting;
    },
    buildCheckpoint: () => {
      checkpointBuilds += 1;
      return { setting, artifact: '<svg id="baseline" />' };
    },
    applyCheckpoint: async (checkpoint) => {
      setting = checkpoint.setting;
    }
  });

  await history.captureBaseline('full baseline');
  await history.runUndoable('Set one', () => {
    setting = 1;
  });
  await history.undo();
  assert.equal(history.getRedoCount(), 1);
  const pending = await history.begin('Pending edit');
  setting = 2;
  const before = history.getDiagnostics();

  await history.initializeIntentBaseline('Loaded session');
  const after = history.getDiagnostics();
  assert.equal(pending.closed, true);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(history.canUndo(), false);
  assert.equal(history.canRedo(), false);
  assert.equal(history.getCurrentCheckpoint(), null);
  assert.equal(history.getCurrentCheckpointSignature(), '');
  assert.equal(after.intentBuilds - before.intentBuilds, 1);
  assert.equal(after.intentSignatureComputations - before.intentSignatureComputations, 1);
  assert.equal(after.artifactCheckpointBuilds - before.artifactCheckpointBuilds, 0);
  assert.equal(
    after.artifactCheckpointSignatureComputations
      - before.artifactCheckpointSignatureComputations,
    0
  );
  assert.equal(after.signatureComputations - before.signatureComputations, 1);
  assert.equal(after.historySvgBytes - before.historySvgBytes, 0);
  assert.equal(checkpointBuilds, 1);

  await history.runUndoable('Set three', () => {
    setting = 3;
  });
  assert.equal(setting, 3);
  await history.undo();
  assert.equal(setting, 2);
  await history.redo();
  assert.equal(setting, 3);
}

{
  let documentState = {
    setting: 'loaded',
    results: [{ name: 'loaded.svg', content: '<svg id="loaded" />' }],
    selectedResultIndex: 0
  };
  let checkpointBuilds = 0;
  const history = createHistoryManager({
    buildIntent: async () => ({ setting: documentState.setting }),
    applyIntent: async (intent) => {
      documentState.setting = intent.setting;
    },
    buildCheckpoint: () => {
      checkpointBuilds += 1;
      return structuredClone(documentState);
    },
    applyCheckpoint: async (checkpoint) => {
      documentState = structuredClone(checkpoint);
    }
  });

  await history.initializeIntentBaseline('Loaded session');
  assert.equal(checkpointBuilds, 0);
  await history.runUndoableCheckpoint('Generate', () => {
    documentState = {
      setting: 'generated',
      results: [
        { name: 'generated-1.svg', content: '<svg id="generated-1" />' },
        { name: 'generated-2.svg', content: '<svg id="generated-2" />' }
      ],
      selectedResultIndex: 1
    };
  });
  assert.equal(checkpointBuilds, 2);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.undoLabel(), 'Generate');
  await history.undo();
  assert.deepEqual(documentState, {
    setting: 'loaded',
    results: [{ name: 'loaded.svg', content: '<svg id="loaded" />' }],
    selectedResultIndex: 0
  });
  await history.redo();
  assert.equal(documentState.results.length, 2);
  assert.equal(documentState.results[1].name, 'generated-2.svg');
  assert.equal(documentState.selectedResultIndex, 1);
}

{
  const captured = [];
  const phases = [];
  const history = createHistoryManager({
    buildIntent: () => ({ setting: 'unchanged' }),
    applyIntent: async () => {},
    buildCheckpoint: () => {
      const checkpoint = { setting: 'unchanged' };
      captured.push(checkpoint);
      return checkpoint;
    },
    applyCheckpoint: async () => {}
  });

  await history.captureBaseline('initial');
  await history.runUndoableCheckpoint('No-op Generate', () => {}, {
    onCheckpointCapture: ({ phase }) => {
      if (phase !== 'before-end' && phase !== 'after-start') return;
      phases.push({ phase, current: history.getCurrentCheckpoint() });
    }
  });

  assert.equal(captured.length, 3);
  assert.deepEqual(phases.map(({ phase }) => phase), ['before-end', 'after-start']);
  assert.equal(phases[0].current, captured[1]);
  assert.equal(phases[1].current, captured[1]);
  assert.equal(history.getCurrentCheckpoint(), captured[2]);
  assert.equal(history.getUndoCount(), 0);
}

{
  let checkpointBuilds = 0;
  const history = createHistoryManager({
    buildIntent: async () => ({ setting: 'loaded' }),
    applyIntent: async () => {},
    buildCheckpoint: () => {
      checkpointBuilds += 1;
      return { results: [{ content: '<svg id="loaded" />' }] };
    },
    applyCheckpoint: async () => {}
  });

  await history.initializeIntentBaseline('Loaded session');
  await assert.rejects(
    history.runUndoableCheckpoint('Failed Generate', async () => {
      throw new Error('generation failed');
    }),
    /generation failed/
  );
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(history.getCurrentCheckpoint(), null);

  const canceled = await history.runUndoableCheckpoint(
    'Canceled Generate',
    async () => ({ status: 'canceled' }),
    { shouldCommit: (result) => result?.status === 'ok' }
  );
  assert.deepEqual(canceled, { status: 'canceled' });
  assert.equal(checkpointBuilds, 2);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(history.getCurrentCheckpoint(), null);
}

{
  let setting = 'loaded-setting';
  let artifact = {
    id: 'loaded',
    identity: { fingerprint: '', compactSignature: 'loaded' },
    retainedBytes: 128,
    fileIds: new Set()
  };
  let workerCalls = 0;
  const restored = [];
  const history = createHistoryManager({
    buildIntent: async () => ({ setting }),
    applyIntent: async (intent) => {
      setting = intent.setting;
    },
    buildCheckpoint: () => {
      throw new Error('Generate replacement must not build a checkpoint.');
    },
    applyCheckpoint: async () => {
      throw new Error('Generate replacement must not apply a checkpoint.');
    },
    captureGeneratedArtifactHandle: () => artifact,
    restoreGeneratedArtifactHandle: async (handle) => {
      artifact = handle;
      restored.push(handle.id);
    },
    compareGeneratedArtifactHandles: (before, after) => Boolean(
      before.identity.fingerprint
      && before.identity.fingerprint === after.identity.fingerprint
      && before.identity.compactSignature === after.identity.compactSignature
    )
  });
  const generated = (id, compactSignature = 'stable') => ({
    id,
    identity: { fingerprint: id.repeat(64).slice(0, 64), compactSignature },
    retainedBytes: 256,
    fileIds: new Set()
  });

  await history.initializeIntentBaseline('Loaded session');
  const beforeFirst = history.getDiagnostics();
  await history.runUndoableArtifactReplacement('Generate A', async (before) => {
    assert.equal(before.id, 'loaded');
    workerCalls += 1;
    artifact = generated('a');
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  const afterFirst = history.getDiagnostics();
  assert.equal(history.getUndoCount(), 1);
  assert.equal(afterFirst.artifactHandleBeforeBuildCount - beforeFirst.artifactHandleBeforeBuildCount, 1);
  assert.equal(afterFirst.artifactHandleAfterBuildCount - beforeFirst.artifactHandleAfterBuildCount, 1);
  assert.equal(afterFirst.artifactFingerprintComparisonCount - beforeFirst.artifactFingerprintComparisonCount, 1);
  assert.equal(afterFirst.artifactReplacementHistoryEntryCount - beforeFirst.artifactReplacementHistoryEntryCount, 1);
  assert.equal(afterFirst.artifactCheckpointBuilds - beforeFirst.artifactCheckpointBuilds, 0);
  assert.equal(afterFirst.artifactCheckpointSignatureComputations - beforeFirst.artifactCheckpointSignatureComputations, 0);
  assert.equal(afterFirst.intentBuilds - beforeFirst.intentBuilds, 0);
  assert.equal(afterFirst.historySvgBytes - beforeFirst.historySvgBytes, 0);
  assert.equal(afterFirst.checkpointEstimatedBytes - beforeFirst.checkpointEstimatedBytes, 0);
  assert.equal(afterFirst.generatedArtifactFullCloneCount - beforeFirst.generatedArtifactFullCloneCount, 0);
  assert.equal(afterFirst.generatedArtifactFullSerializationCount - beforeFirst.generatedArtifactFullSerializationCount, 0);
  assert.equal(
    afterFirst.manualCancelFullArtifactSnapshotBuildCount
      - beforeFirst.manualCancelFullArtifactSnapshotBuildCount,
    0
  );

  await history.undo();
  assert.equal(artifact.id, 'loaded');
  await history.redo();
  assert.equal(artifact.id, 'a');
  assert.equal(workerCalls, 1, 'Undo/Redo must not invoke generation');

  const beforeUnchanged = history.getDiagnostics();
  const undoCountBeforeUnchanged = history.getUndoCount();
  const redoCountBeforeUnchanged = history.getRedoCount();
  await history.runUndoableArtifactReplacement('Generate A unchanged', async () => {
    workerCalls += 1;
    artifact = generated('a');
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  const afterUnchanged = history.getDiagnostics();
  assert.equal(history.getUndoCount(), undoCountBeforeUnchanged);
  assert.equal(history.getRedoCount(), redoCountBeforeUnchanged);
  assert.equal(
    afterUnchanged.artifactReplacementHistoryEntryCount
      - beforeUnchanged.artifactReplacementHistoryEntryCount,
    0
  );
  assert.equal(afterUnchanged.artifactHandleBeforeBuildCount - beforeUnchanged.artifactHandleBeforeBuildCount, 1);
  assert.equal(afterUnchanged.artifactHandleAfterBuildCount - beforeUnchanged.artifactHandleAfterBuildCount, 1);
  assert.equal(
    afterUnchanged.artifactFingerprintComparisonCount
      - beforeUnchanged.artifactFingerprintComparisonCount,
    1
  );

  const pendingSetting = await history.begin('Change render setting');
  setting = 'changed-setting';
  await history.runUndoableArtifactReplacement('Generate B', async () => {
    workerCalls += 1;
    artifact = generated('b');
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  assert.equal(pendingSetting.closed, true);
  assert.equal(pendingSetting.deferAdapterCommit, true);
  assert.equal(history.undoLabel(), 'Generate B');
  await history.undo();
  assert.equal(artifact.id, 'a');
  assert.equal(setting, 'changed-setting');
  assert.equal(history.undoLabel(), 'Change render setting');
  await history.undo();
  assert.equal(setting, 'loaded-setting');
  await history.redo();
  await history.redo();
  assert.equal(setting, 'changed-setting');
  assert.equal(artifact.id, 'b');

  await history.runUndoableArtifactReplacement('Generate C', async () => {
    workerCalls += 1;
    artifact = generated('c');
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  await history.undo();
  assert.equal(artifact.id, 'b');
  await history.undo();
  assert.equal(artifact.id, 'a');
  await history.redo();
  assert.equal(artifact.id, 'b');
  await history.redo();
  assert.equal(artifact.id, 'c');

  const countBeforeFailure = history.getUndoCount();
  const failed = await history.runUndoableArtifactReplacement('Failed Generate', async (before) => {
    artifact = generated('d');
    artifact = before;
    return { status: 'error' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  assert.deepEqual(failed, { status: 'error' });
  assert.equal(artifact.id, 'c');
  assert.equal(history.getUndoCount(), countBeforeFailure);
  assert.equal(history.getRedoCount(), 0);

  await assert.rejects(
    history.runUndoableArtifactReplacement('Thrown Generate', async () => {
      artifact = generated('d');
      throw new Error('injected generation failure');
    }),
    /injected generation failure/
  );
  assert.equal(artifact.id, 'c');
  assert.equal(restored.at(-1), 'c');
}

{
  let artifact = {
    identity: { fingerprint: 'a'.repeat(64), compactSignature: 'a' },
    retainedBytes: 80,
    fileIds: []
  };
  const history = createHistoryManager({
    maxBytes: 100,
    buildIntent: async () => ({}),
    applyIntent: async () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: async () => {},
    captureGeneratedArtifactHandle: () => artifact,
    restoreGeneratedArtifactHandle: async (handle) => {
      artifact = handle;
    }
  });
  await history.initializeIntentBaseline('Loaded session');
  await history.runUndoableArtifactReplacement('Oversized Generate', async () => {
    artifact = {
      identity: { fingerprint: 'b'.repeat(64), compactSignature: 'b' },
      retainedBytes: 80,
      fileIds: []
    };
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  assert.equal(history.getDiagnostics().artifactReplacementHistoryEntryCount, 1);
  assert.equal(history.canUndo(), false);
  assert.match(history.historyLimitMessage.value, /discarded/);
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
  const circularPrimary = makeFile('circular-primary.gb', 101);
  const circularDepth = makeFile('circular-depth.tsv', 102);
  const conservationBlast = makeFile('derived-conservation.tsv', 103);
  const conservationFasta = makeFile('conservation-subject.fa', 104);
  const conservationSequenceSource = makeFile('derived-conservation-source.fa', 105);
  const linearPrimary = makeFile('linear-primary.gb', 106);
  const linearDepth = makeFile('linear-depth.tsv', 107);
  const linearComparison = makeFile('linear-comparison.tsv', 108);
  const canonicalProtein = makeFile('canonical-protein.tsv', 109);
  const canonicalOrthogroups = makeFile('canonical-orthogroups.json', 110);
  const canonicalCollinearity = makeFile('canonical-collinearity.json', 111);
  const unreferenced = makeFile('unreferenced.tmp', 112);
  const retainedFiles = [
    circularPrimary,
    circularDepth,
    conservationBlast,
    conservationFasta,
    conservationSequenceSource,
    linearPrimary,
    linearDepth,
    linearComparison,
    canonicalProtein,
    canonicalOrthogroups,
    canonicalCollinearity
  ];
  const retainedIds = retainedFiles.map((file) => fileStore.describeFile(file).fileId);
  const unreferencedId = fileStore.describeFile(unreferenced).fileId;
  const state = {
    files: {
      c_gb: circularPrimary,
      c_gff: null,
      c_fasta: null,
      c_depth: [circularDepth],
      c_conservation_blasts: [conservationBlast],
      c_conservation_blasts_source: 'losat-cache',
      c_conservation_fastas: [conservationFasta],
      c_conservation_sequence_sources: [conservationSequenceSource],
      linearCanonicalComparisons: [
        { kind: 'precomputedProteinComparison', file: canonicalProtein },
        { kind: 'orthogroupResult', file: canonicalOrthogroups },
        { kind: 'collinearityResult', file: canonicalCollinearity }
      ],
      d_color: null,
      t_color: null,
      blacklist: null,
      whitelist: null,
      qualifier_priority: null
    },
    linearSeqs: [{
      uid: 'linear-a',
      gb: linearPrimary,
      gff: null,
      fasta: null,
      depth: [linearDepth]
    }],
    linearComparisonPlan: {
      mode: 'selected',
      defaultSource: 'upload',
      edges: [{ id: 'edge-a-b', file: linearComparison }]
    }
  };
  let setting = 'loaded';
  let artifactBuilds = 0;
  const snapshots = createHistorySnapshotService({ state, fileStore });
  const history = createHistoryManager({
    fileStore,
    collectCurrentFileIds: snapshots.collectCurrentFileIds,
    buildIntent: async () => ({ setting }),
    applyIntent: async (intent) => {
      setting = intent.setting;
    },
    buildCheckpoint: () => {
      artifactBuilds += 1;
      throw new Error('lightweight baseline must not build an artifact checkpoint');
    },
    applyCheckpoint: snapshots.applyArtifactCheckpoint
  });

  await history.initializeIntentBaseline('Loaded session');
  assert.equal(artifactBuilds, 0);
  retainedIds.forEach((fileId) => assert.equal(fileStore.has(fileId), true));
  assert.equal(fileStore.has(unreferencedId), false);

  await history.runUndoable('Change ordinary setting', () => {
    setting = 'edited';
  });
  await history.undo();
  assert.equal(setting, 'loaded');
  retainedIds.forEach((fileId) => assert.equal(fileStore.has(fileId), true));
  await history.redo();
  assert.equal(setting, 'edited');
  retainedIds.forEach((fileId) => assert.equal(fileStore.has(fileId), true));
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

{
  const fileStore = createHistoryFileStore();
  const resultA = { name: 'a.svg', content: '<svg id="a" />' };
  const extractedA = [{ svg_id: 'feature-a', payload: 'x'.repeat(100_000) }];
  const biologicalA = [{ biological_feature_id: 'bio-a', payload: 'y'.repeat(100_000) }];
  const catalogA = { schema: 3, items: [{ payload: 'z'.repeat(100_000) }] };
  const groupsA = [{ id: 'group-a', members: ['feature-a'] }];
  const proteinIdentityManifestA = { schema: 2, records: [{ id: 'manifest-a' }] };
  const legacyProteinRawCandidatesA = { schema: 1, entries: [{ id: 'raw-a' }] };
  const legacyProteinDerivedEvidenceA = {
    schema: 1,
    entries: [{ id: 'derived-evidence-a' }]
  };
  const rawCacheValueA = { text: 'raw-a\n' };
  const derivedCacheValueA = { payload: { id: 'derived-a' } };
  const losatCacheInfoA = [{ key: 'raw-a' }];
  const appliedPaletteColorsA = { CDS: '#54bcf8' };
  const pendingPaletteColorsA = {};
  const state = {
    files: {},
    linearSeqs: [],
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] },
    mode: ref('circular'),
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    results: ref([resultA]),
    selectedResultIndex: ref(0),
    extractedFeatures: ref(extractedA),
    biologicalFeatures: ref(biologicalA),
    featureCatalog: ref(catalogA),
    featureSelectorSafetyScope: ref([]),
    featureRecordIds: ref(['record-a']),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: { 'feature-a': '#112233' },
    featureVisibilityManualRules: [],
    featureVisibilityOverrides: {},
    featureVisibilitySelectorCache: {
      'feature-a': { qualifier: 'hash', value: 'feature-a' }
    },
    labelTextFeatureOverrides: {},
    canonicalLabelOverrideRows: ref([]),
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    labelOverrideContextKey: ref('context-a'),
    legendEntries: ref([{ caption: 'A', color: '#112233' }]),
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref(['A']),
    originalLegendColors: ref({ A: '#112233' }),
    legendColorOverrides: {},
    legendStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()),
    featureStrokeOverrides: {},
    originalSvgStroke: ref({ color: '#000000', width: 1 }),
    orthogroups: ref(groupsA),
    featureOrthogroupIndex: ref(new Map([['feature-a', 'group-a']])),
    selectedOrthogroupId: ref('group-a'),
    selectedOrthogroupAlignmentFeature: ref('feature-a'),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    collinearGroups: ref([]),
    lastRunInfo: ref({ resultCount: 1 }),
    pairwiseMatchFactors: ref({}),
    trackSlotResolvedGeometry: ref({ schema: 1 }),
    appliedPaletteName: ref('default'),
    appliedPaletteColors: ref(appliedPaletteColorsA),
    pendingPaletteName: ref(''),
    pendingPaletteColors: ref(pendingPaletteColorsA),
    resultGenerationKey: ref(1),
    resultPanelTab: ref('preview'),
    zoom: ref(1),
    canvasPan: { x: 0, y: 0 },
    errorLog: ref(null),
    editableLabels: ref([]),
    labelTextScopeDialog: {},
    featureEditorStatus: { status: 'summary-ready', summaryCount: 1 },
    featureExtractionPending: ref(false),
    featureExtractionError: ref(null),
    labelOverrideBuildWarning: ref(''),
    proteinIdentityManifest: ref(proteinIdentityManifestA),
    legacyProteinRawCandidates: ref(legacyProteinRawCandidatesA),
    legacyProteinDerivedEvidence: ref(legacyProteinDerivedEvidenceA),
    losatCache: ref(new Map([['raw-a', rawCacheValueA]])),
    losatDerivedCache: ref(new Map([['derived-a', derivedCacheValueA]])),
    losatCacheInfo: ref(losatCacheInfoA),
    matchSequenceRegistry: (() => {
      let owner = Object.freeze({
        sources: new Map([['record-a', {
          key: 'record-a', recordId: 'record-a', aliases: [], sequence: 'ACGT'
        }]]),
        ambiguousKeys: new Set()
      });
      return {
        captureTrustedOwner: () => owner,
        replaceTrustedOwner: (nextOwner) => { owner = nextOwner; },
        values: () => Array.from(owner.sources.values())
      };
    })(),
    skipCaptureBaseConfig: ref(false),
    skipPositionReapply: ref(false),
    skipExtractOnSvgChange: ref(false),
    trustedArtifactRestoreInProgress: ref(false),
    semanticFileWatchersSuppressed: ref(false)
  };
  const captureAuthorityOwners = () => ({
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
    legacyProteinDerivedEvidence: state.legacyProteinDerivedEvidence.value,
    losatCache: state.losatCache.value,
    losatDerivedCache: state.losatDerivedCache.value,
    losatCacheInfo: state.losatCacheInfo.value
  });
  const assertAuthorityOwners = (expected) => {
    assert.equal(state.proteinIdentityManifest.value, expected.proteinIdentityManifest);
    assert.equal(state.legacyProteinRawCandidates.value, expected.legacyProteinRawCandidates);
    assert.equal(
      state.legacyProteinDerivedEvidence.value,
      expected.legacyProteinDerivedEvidence
    );
    assert.equal(state.losatCache.value, expected.losatCache);
    assert.equal(state.losatDerivedCache.value, expected.losatDerivedCache);
    assert.equal(state.losatCacheInfo.value, expected.losatCacheInfo);
  };
  const authorityOwnersA = captureAuthorityOwners();
  let resultAdmissionCalls = 0;
  const snapshots = createHistorySnapshotService({
    state,
    fileStore,
    nextTick: async () => {},
    buildUiStateData: () => ({
      mode: state.mode.value,
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      selectedResultIndex: state.selectedResultIndex.value,
      zoom: state.zoom.value,
      canvasPan: { ...state.canvasPan },
      appliedPaletteName: state.appliedPaletteName.value,
      appliedPaletteColors: { ...state.appliedPaletteColors.value },
      pendingPaletteName: state.pendingPaletteName.value,
      pendingPaletteColors: { ...state.pendingPaletteColors.value }
    }),
    applyUiStateData: (ui) => {
      state.mode.value = ui.mode;
      state.selectedResultIndex.value = ui.selectedResultIndex;
    },
    serializeResults: () => {
      throw new Error('Generated Artifact Handles must not serialize Results.');
    },
    applyResultsData: () => {
      resultAdmissionCalls += 1;
      throw new Error('Trusted artifact restore must not use Result admission.');
    },
    applyFeatureStateData: (features) => {
      state.extractedFeatures.value = features.extractedFeatures;
      state.biologicalFeatures.value = features.biologicalFeatures;
      state.featureSelectorSafetyScope.value = features.featureSelectorSafetyScope;
      state.featureRecordIds.value = features.featureRecordIds;
      state.selectedFeatureRecordIdx.value = features.selectedFeatureRecordIdx;
      Object.assign(state.featureColorOverrides, features.featureColorOverrides);
    },
    applyEditorStateData: (editorState, options) => {
      assert.equal(options.normalized, true);
      assert.equal(options.adoptCatalog, true);
      state.legendEntries.value = editorState.legend.entries;
      state.featureCatalog.value = editorState.featureCatalog;
    },
    buildRunStateData: () => ({
      lastRunInfo: structuredClone(state.lastRunInfo.value),
      pairwiseMatchFactors: structuredClone(state.pairwiseMatchFactors.value)
    }),
    applyRunStateData: (runState) => {
      state.lastRunInfo.value = structuredClone(runState.lastRunInfo);
      state.pairwiseMatchFactors.value = structuredClone(runState.pairwiseMatchFactors);
    }
  });
  snapshots.setGeneratedArtifactIdentity({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint: 'a'.repeat(64),
    retainedBytes: 400_000
  }, { results: state.results.value });

  const handleA = snapshots.captureGeneratedArtifactHandle();
  assert.equal(Object.isFrozen(handleA), true);
  assert.equal(handleA.ownerSet.results, state.results.value);
  assert.equal(handleA.ownerSet.results[0], resultA);
  assert.equal(handleA.ownerSet.extractedFeatures, extractedA);
  assert.equal(handleA.ownerSet.biologicalFeatures, biologicalA);
  assert.equal(handleA.ownerSet.featureCatalog, catalogA);
  assert.equal(handleA.ownerSet.orthogroups, groupsA);
  assert.equal(handleA.ownerSet.proteinIdentityManifest, proteinIdentityManifestA);
  assert.equal(handleA.ownerSet.legacyProteinRawCandidates, legacyProteinRawCandidatesA);
  assert.equal(handleA.ownerSet.legacyProteinDerivedEvidence, legacyProteinDerivedEvidenceA);
  assert.equal(handleA.ownerSet.losatCache, authorityOwnersA.losatCache);
  assert.equal(handleA.ownerSet.losatDerivedCache, authorityOwnersA.losatDerivedCache);
  assert.equal(handleA.ownerSet.losatCacheInfo, losatCacheInfoA);
  assert.equal(handleA.ownerSet.appliedPaletteColors, appliedPaletteColorsA);
  assert.equal(handleA.ownerSet.pendingPaletteColors, pendingPaletteColorsA);
  assert.equal(
    Object.fromEntries(handleA.mutableIntent.features.featureVisibilitySelectorCache)['feature-a'].value,
    'feature-a'
  );
  assert.ok(handleA.retainedBytes >= 400_000);

  state.resultGenerationKey.value += 1;
  state.zoom.value = 2;
  state.canvasPan.x = 120;
  const transientOnlyHandle = snapshots.captureGeneratedArtifactHandle();
  assert.equal(
    snapshots.compareGeneratedArtifactHandles(handleA, transientOnlyHandle),
    true,
    'Generation counters and preview viewport state must not create artifact History.'
  );
  state.resultGenerationKey.value = 1;
  state.zoom.value = 1;
  state.canvasPan.x = 0;

  const originalOwnerSet = snapshots.captureGeneratedArtifactOwnerSet();
  const poisonResult = { name: 'poison.svg' };
  Object.defineProperty(poisonResult, 'content', {
    get() { throw new Error('Result content must not be read during capture.'); }
  });
  const poisonCatalog = {};
  Object.defineProperty(poisonCatalog, 'items', {
    get() { throw new Error('Catalog rows must not be read during capture.'); }
  });
  const poisonCache = new Proxy(new Map(), {
    get() { throw new Error('Cache properties and iterators must not be read during capture.'); }
  });
  const poisonSequenceOwner = {};
  Object.defineProperty(poisonSequenceOwner, 'sources', {
    get() { throw new Error('Sequence owners must not be traversed during capture.'); }
  });
  const poisonResource = { name: 'poison.gb', size: 123, type: 'text/plain' };
  Object.defineProperty(poisonResource, 'data', {
    get() { throw new Error('Resource payloads must not be read during capture.'); }
  });
  const structuralMetrics = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => structuralMetrics.push(metric)
  };
  state.results.value = [poisonResult];
  state.featureCatalog.value = poisonCatalog;
  state.losatCache.value = poisonCache;
  state.losatDerivedCache.value = poisonCache;
  state.matchSequenceRegistry.replaceTrustedOwner(poisonSequenceOwner);
  state.files.c_gb = poisonResource;
  snapshots.setGeneratedArtifactRuntimeOwner({
    capture: () => ({ retainedBytes: 999, payload: poisonResource }),
    restore: async () => {}
  });
  snapshots.setGeneratedArtifactIdentity({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint: 'c'.repeat(64),
    retainedBytes: 123_456,
    resultBytes: 654_321
  }, { results: state.results.value });
  const poisonHandle = snapshots.captureGeneratedArtifactHandle();
  assert.equal(poisonHandle.ownerSet.results[0], poisonResult);
  assert.equal(poisonHandle.ownerSet.featureCatalog, poisonCatalog);
  assert.equal(poisonHandle.ownerSet.losatCache, poisonCache);
  assert.equal(poisonHandle.ownerSet.matchSequenceOwner, poisonSequenceOwner);
  assert.equal(
    poisonHandle.retainedBytes,
    123_456 + 654_321 + 999 + poisonHandle.identity.compactSignature.length * 2
  );
  assert.deepEqual(
    structuralMetrics.filter(({ name }) => name === 'generatedArtifactFullSignatureCount')
      .map(({ value }) => value),
    [0]
  );
  assert.deepEqual(
    structuralMetrics.filter(({ name }) => name === 'generatedArtifactHeavyTraversalCount')
      .map(({ value }) => value),
    [0]
  );
  assert.deepEqual(
    structuralMetrics.filter(({ name }) => name === 'generatedArtifactMutableIntentSnapshotCount')
      .map(({ value }) => value),
    [1]
  );
  let installedResults = null;
  snapshots.installGeneratedArtifactOwnerSet(originalOwnerSet, {
    installResults: (nextResults) => {
      installedResults = nextResults;
      state.results.value = nextResults;
    }
  });
  assert.equal(installedResults, originalOwnerSet.results);
  assert.equal(state.results.value, originalOwnerSet.results);
  state.files.c_gb = null;
  snapshots.setGeneratedArtifactRuntimeOwner({
    capture: () => null,
    restore: async () => {}
  });
  snapshots.setGeneratedArtifactIdentity({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint: 'a'.repeat(64),
    retainedBytes: 400_000
  }, { results: state.results.value });
  globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;

  state.results.value = [{ name: 'a.svg', content: '<svg id="b" />' }];
  const invalidatedIdentityHandle = snapshots.captureGeneratedArtifactHandle();
  assert.equal(
    invalidatedIdentityHandle.identity.fingerprint,
    '',
    'Replacing a Result object must invalidate its Worker-byte identity even at equal length.'
  );
  state.results.value = [{ name: 'edited.svg', content: '<svg id="edited" />' }];
  state.featureColorOverrides['feature-a'] = '#abcdef';
  state.featureVisibilitySelectorCache['feature-a'] = {
    qualifier: 'hash',
    value: 'edited-feature'
  };
  state.legendEntries.value[0].caption = 'Edited current legend';
  state.extractedFeatures.value = [{ svg_id: 'feature-b' }];
  state.biologicalFeatures.value = [{ biological_feature_id: 'bio-b' }];
  state.featureCatalog.value = { schema: 3, items: [] };
  state.orthogroups.value = [{ id: 'group-b' }];
  const proteinIdentityManifestB = { schema: 2, records: [{ id: 'manifest-b' }] };
  const legacyProteinRawCandidatesB = { schema: 1, entries: [{ id: 'raw-b' }] };
  const legacyProteinDerivedEvidenceB = {
    schema: 1,
    entries: [{ id: 'derived-evidence-b' }]
  };
  const rawCacheValueB = { text: 'raw-b\n' };
  const derivedCacheValueB = { payload: { id: 'derived-b' } };
  const losatCacheInfoB = [{ key: 'raw-b' }];
  state.proteinIdentityManifest.value = proteinIdentityManifestB;
  state.legacyProteinRawCandidates.value = legacyProteinRawCandidatesB;
  state.legacyProteinDerivedEvidence.value = legacyProteinDerivedEvidenceB;
  state.losatCache.value = new Map([['raw-b', rawCacheValueB]]);
  state.losatDerivedCache.value = new Map([['derived-b', derivedCacheValueB]]);
  state.losatCacheInfo.value = losatCacheInfoB;
  const appliedPaletteColorsB = { ...appliedPaletteColorsA };
  const pendingPaletteColorsB = { ...pendingPaletteColorsA };
  state.appliedPaletteColors.value = appliedPaletteColorsB;
  state.pendingPaletteColors.value = pendingPaletteColorsB;
  const authorityOwnersB = captureAuthorityOwners();
  snapshots.setGeneratedArtifactIdentity({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint: 'b'.repeat(64),
    retainedBytes: 1024
  }, { results: state.results.value });
  const handleB = snapshots.captureGeneratedArtifactHandle();

  assert.equal(handleA.ownerSet.results[0], resultA);
  assert.equal(handleA.mutableIntent.features.featureColorOverrides['feature-a'], '#112233');
  assert.equal(handleA.mutableIntent.editorState.legend.entries[0].caption, 'A');
  assert.equal(
    Object.fromEntries(handleA.mutableIntent.features.featureVisibilitySelectorCache)['feature-a'].value,
    'feature-a'
  );
  assert.equal(handleB.ownerSet.results[0].name, 'edited.svg');
  assert.equal(handleB.ownerSet.proteinIdentityManifest, proteinIdentityManifestB);
  assert.equal(handleB.ownerSet.legacyProteinRawCandidates, legacyProteinRawCandidatesB);
  assert.equal(handleB.ownerSet.legacyProteinDerivedEvidence, legacyProteinDerivedEvidenceB);
  assert.equal(handleB.ownerSet.losatCache, authorityOwnersB.losatCache);
  assert.equal(handleB.ownerSet.losatDerivedCache, authorityOwnersB.losatDerivedCache);
  assert.equal(handleB.ownerSet.losatCacheInfo, losatCacheInfoB);
  assert.equal(handleB.ownerSet.appliedPaletteColors, appliedPaletteColorsB);
  assert.equal(handleB.ownerSet.pendingPaletteColors, pendingPaletteColorsB);

  state.proteinIdentityManifest.value = { schema: 2, records: [{ id: 'replacement' }] };
  state.legacyProteinRawCandidates.value = { schema: 1, entries: [] };
  state.legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
  state.losatCache.value = new Map([['replacement', { text: 'replacement\n' }]]);
  state.losatDerivedCache.value = new Map([['replacement', { payload: {} }]]);
  state.losatCacheInfo.value = [{ key: 'replacement' }];
  await snapshots.restoreGeneratedArtifactHandle(handleA);
  assert.equal(resultAdmissionCalls, 0);
  assert.equal(state.results.value[0], resultA);
  assert.equal(state.extractedFeatures.value, extractedA);
  assert.equal(state.biologicalFeatures.value, biologicalA);
  assert.equal(state.featureCatalog.value, catalogA);
  assert.equal(state.orthogroups.value, groupsA);
  assert.equal(state.featureColorOverrides['feature-a'], '#112233');
  assert.equal(state.legendEntries.value[0].caption, 'A');
  assert.equal(state.featureVisibilitySelectorCache['feature-a'].value, 'feature-a');
  assertAuthorityOwners(authorityOwnersA);
  assert.equal(state.appliedPaletteColors.value, appliedPaletteColorsA);
  assert.equal(state.pendingPaletteColors.value, pendingPaletteColorsA);
  const recapturedA = snapshots.captureGeneratedArtifactHandle();
  assert.equal(recapturedA.identity.fingerprint, 'a'.repeat(64));
  assert.equal(recapturedA.retainedBytes, handleA.retainedBytes);

  state.results.value = [{ ...state.results.value[0], content: '<svg id="normal-edit" />' }];
  state.featureColorOverrides['feature-a'] = '#ffffff';
  state.legendEntries.value[0].caption = 'Normal supported edit';
  assert.equal(handleA.ownerSet.results[0].content, '<svg id="a" />');
  assert.equal(handleA.mutableIntent.features.featureColorOverrides['feature-a'], '#112233');
  assert.equal(handleA.mutableIntent.editorState.legend.entries[0].caption, 'A');

  await snapshots.restoreGeneratedArtifactHandle(handleB);
  assert.equal(state.results.value[0].name, 'edited.svg');
  assert.equal(state.extractedFeatures.value[0].svg_id, 'feature-b');
  assertAuthorityOwners(authorityOwnersB);
  assert.equal(state.appliedPaletteColors.value, appliedPaletteColorsB);
  assert.equal(state.pendingPaletteColors.value, pendingPaletteColorsB);

  await Promise.all([
    snapshots.restoreGeneratedArtifactHandle(handleA),
    snapshots.restoreGeneratedArtifactHandle(handleB)
  ]);
  assert.equal(state.trustedArtifactRestoreInProgress.value, false);
  assert.equal(state.semanticFileWatchersSuppressed.value, false);

  await snapshots.restoreGeneratedArtifactHandle(handleA);
  const installAuthorityOwnersB = () => {
    snapshots.installGeneratedArtifactOwnerSet(handleB.ownerSet, {
      selectedResultIndex: handleB.mutableIntent.ui.selectedResultIndex
    });
    snapshots.setGeneratedArtifactIdentity({
      schema: 1,
      algorithm: 'SHA-256',
      fingerprint: 'b'.repeat(64),
      retainedBytes: 1024
    }, { results: state.results.value });
  };
  const ownerHistory = createHistoryManager({
    fileStore,
    buildIntent: async () => ({}),
    applyIntent: async () => {},
    buildCheckpoint: () => {
      throw new Error('LOSAT authority History must not build a full checkpoint.');
    },
    applyCheckpoint: async () => {
      throw new Error('LOSAT authority History must not apply a full checkpoint.');
    },
    captureGeneratedArtifactHandle: snapshots.captureGeneratedArtifactHandle,
    restoreGeneratedArtifactHandle: snapshots.restoreGeneratedArtifactHandle,
    compareGeneratedArtifactHandles: snapshots.compareGeneratedArtifactHandles
  });
  await ownerHistory.initializeIntentBaseline('LOSAT authority baseline');
  const replacementLifecycle = [];
  const replacementMetrics = [];
  const hooksBeforeReplacement = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onSessionLifecycleEvent: ({ name }) => replacementLifecycle.push(name),
    onStructuralMetric: (metric) => replacementMetrics.push(metric)
  };
  let successfulGenerateCalls = 0;
  await ownerHistory.runUndoableArtifactReplacement('Generate LOSAT owners B', async () => {
    successfulGenerateCalls += 1;
    installAuthorityOwnersB();
    return { status: 'ok' };
  }, { shouldCommit: (result) => result.status === 'ok' });
  assertAuthorityOwners(authorityOwnersB);
  assert.equal(ownerHistory.getUndoCount(), 1);
  assert.deepEqual(replacementLifecycle.slice(0, 2), [
    'history.before-capture-started',
    'history.before-capture-completed'
  ]);
  assert.equal(
    replacementMetrics.filter(({ name }) => name === 'historyReplacementCount')
      .reduce((total, { value }) => total + value, 0),
    1
  );
  await ownerHistory.runUndoableArtifactReplacement(
    'No-op Generate',
    async () => ({ status: 'ok' }),
    { shouldCommit: (result) => result.status === 'ok' }
  );
  assert.equal(ownerHistory.getUndoCount(), 1);
  assert.equal(
    replacementMetrics.filter(({ name }) => name === 'historyReplacementCount')
      .reduce((total, { value }) => total + value, 0),
    1
  );
  await ownerHistory.undo();
  assertAuthorityOwners(authorityOwnersA);
  await ownerHistory.redo();
  assertAuthorityOwners(authorityOwnersB);
  assert.equal(successfulGenerateCalls, 1, 'Undo/Redo must not rerun Generate.');

  const restoredHandleBeforePaletteReplacement = snapshots.captureGeneratedArtifactHandle();
  const equivalentAppliedPaletteColors = { ...state.appliedPaletteColors.value };
  const equivalentPendingPaletteColors = { ...state.pendingPaletteColors.value };
  state.appliedPaletteColors.value = equivalentAppliedPaletteColors;
  state.pendingPaletteColors.value = equivalentPendingPaletteColors;
  const equivalentPaletteHandle = snapshots.captureGeneratedArtifactHandle();
  assert.notEqual(equivalentAppliedPaletteColors, handleB.ownerSet.appliedPaletteColors);
  assert.notEqual(equivalentPendingPaletteColors, handleB.ownerSet.pendingPaletteColors);
  assert.equal(equivalentPaletteHandle.identity.fingerprint, 'b'.repeat(64));
  assert.equal(
    snapshots.compareGeneratedArtifactHandles(
      restoredHandleBeforePaletteReplacement,
      equivalentPaletteHandle
    ),
    true,
    'Equivalent palette owner replacements must not change semantic artifact identity.'
  );
  await ownerHistory.runUndoableArtifactReplacement(
    'No-op Generate after Undo/Redo',
    async () => ({ status: 'ok' }),
    { shouldCommit: (result) => result.status === 'ok' }
  );
  assert.equal(ownerHistory.getUndoCount(), 1);
  assert.equal(
    replacementMetrics.filter(({ name }) => name === 'historyReplacementCount')
      .reduce((total, { value }) => total + value, 0),
    1
  );

  const changedAppliedPaletteColors = {
    ...equivalentAppliedPaletteColors,
    CDS: '#abcdef'
  };
  await ownerHistory.runUndoableArtifactReplacement(
    'Change applied palette intent',
    async () => {
      state.appliedPaletteName.value = 'custom';
      state.appliedPaletteColors.value = changedAppliedPaletteColors;
      return { status: 'ok' };
    },
    { shouldCommit: (result) => result.status === 'ok' }
  );
  assert.equal(ownerHistory.getUndoCount(), 2);
  assert.equal(
    replacementMetrics.filter(({ name }) => name === 'historyReplacementCount')
      .reduce((total, { value }) => total + value, 0),
    2,
    'A real bounded palette intent change must create History replacement.'
  );
  await ownerHistory.undo();
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.equal(state.appliedPaletteColors.value, equivalentAppliedPaletteColors);
  await ownerHistory.redo();
  assert.equal(state.appliedPaletteName.value, 'custom');
  assert.equal(state.appliedPaletteColors.value, changedAppliedPaletteColors);
  globalThis.__GBDRAW_TEST_HOOKS__ = hooksBeforeReplacement;

  const installDiscardedAuthorityOwners = (status) => {
    state.proteinIdentityManifest.value = { schema: 2, records: [{ id: status }] };
    state.legacyProteinRawCandidates.value = { schema: 1, entries: [{ id: status }] };
    state.legacyProteinDerivedEvidence.value = { schema: 1, entries: [{ id: status }] };
    state.losatCache.value = new Map([[status, { text: `${status}\n` }]]);
    state.losatDerivedCache.value = new Map([[status, { payload: { id: status } }]]);
    state.losatCacheInfo.value = [{ key: status }];
  };
  const undoCountBeforeDiscardedRuns = ownerHistory.getUndoCount();
  await assert.rejects(
    ownerHistory.runUndoableArtifactReplacement('Failed LOSAT Generate', async () => {
      installDiscardedAuthorityOwners('failed');
      throw new Error('injected LOSAT generation failure');
    }),
    /injected LOSAT generation failure/
  );
  assertAuthorityOwners(authorityOwnersB);

  for (const status of ['error', 'canceled', 'stale']) {
    const result = await ownerHistory.runUndoableArtifactReplacement(
      `${status} LOSAT Generate`,
      async (before) => {
        installDiscardedAuthorityOwners(status);
        await snapshots.restoreGeneratedArtifactHandle(before);
        return { status };
      },
      { shouldCommit: (outcome) => outcome.status === 'ok' }
    );
    assert.deepEqual(result, { status });
    assertAuthorityOwners(authorityOwnersB);
  }
  assert.equal(ownerHistory.getUndoCount(), undoCountBeforeDiscardedRuns);
}

console.log('history tests passed');
