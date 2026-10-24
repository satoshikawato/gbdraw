import assert from 'node:assert/strict';
import { cp, mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import vm from 'node:vm';
import { withDrawings, withModeDrawings } from './helpers/drawing-state.mjs';

// The fixture as both drawings (PD-OI-090: one History per drawing), with a
// Linear view whose own `mode` is 'linear', so each drawing's intent holds the
// bindings of its own slots; `state.mode` picks the shown drawing.
const withTwoDrawings = (state) => withModeDrawings(state, {
  circular: state, linear: Object.create(state, { mode: { value: 'linear' } })
});

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'services');
const webSourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-history-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'web-ux-profile.js'),
  await readFile(join(webSourceDir, 'web-ux-profile.js'), 'utf8'),
  'utf8'
);
await cp(webSourceDir, tempDir, {recursive:true});
await mkdir(join(tempDir, 'services'), { recursive: true });
await mkdir(join(tempDir, 'app'), { recursive: true });
for (const filename of [
  'canonical-comparisons.js',
  'history.js',
  'history-files.js',
  'history-snapshot.js',
  'json-clone.js',
  'feature-selector.js',
  'feature-utils.js',
  'feature-visibility.js',
  'layout-preferences.js',
  'runtime-test-hooks.js',
  'svg-serialization.js'
]) {
  await writeFile(
    join(tempDir, 'services', filename),
    await readFile(join(sourceDir, filename), 'utf8'),
    'utf8'
  );
}

const { createHistoryManager, createHistoryRouter } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));
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

  // E1: an Undo or Redo is pending from its call, before it restores, to its
  // end; the diagram mode switch waits for it.
  assert.equal(history.traversalPending(), false);
  const undoing = history.undo();
  assert.equal(history.traversalPending(), true);
  assert.equal(history.restoring.value, false);
  await undoing;
  assert.equal(history.traversalPending(), false);
  const redoing = history.redo();
  assert.equal(history.traversalPending(), true);
  await redoing;
  assert.equal(history.traversalPending(), false);
}

// E1 (REVIEW-2 P-a): the composition root may report one Undo or Redo step
// busy (a step that switches the diagram mode while a label rerender runs);
// that step is not applied and stays on its stack.
{
  let value = 0;
  /** @type {{ status: string, reason: string } | null} */
  let busy = null;
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => { value = snapshot.value; };
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState,
    stepAvailability: (step) => (step.type === 'intent' ? busy : null)
  });
  await history.captureBaseline('initial');
  await history.runUndoable('Increment', () => { value = 1; });
  busy = { status: 'busy', reason: 'Updating diagram. Retry after the update finishes.' };
  assert.deepEqual(await history.undo(), busy);
  assert.equal(value, 1);
  assert.equal(history.getUndoCount(), 1);
  busy = null;
  assert.equal(await history.undo(), true);
  assert.equal(value, 0);
  busy = { status: 'busy', reason: 'Updating diagram. Retry after the update finishes.' };
  assert.deepEqual(await history.redo(), busy);
  assert.equal(value, 0);
  assert.equal(history.getRedoCount(), 1);
}

// OV-1339: the enabled state reads the availability the action reads. While
// the top step is busy, Undo (or Redo) and its label are off; they come back,
// reactively and without another History change, when the busy state ends.
// The header buttons' computed is built as app-setup.js builds it, over the
// vendored Vue.
{
  const vueContext = vm.createContext({});
  vm.runInContext(
    await readFile(join(repoRoot, 'gbdraw', 'web', 'vendor', 'vue', 'vue.global.prod.js'), 'utf8'),
    vueContext
  );
  const Vue = vm.runInContext('Vue', vueContext);
  let value = 0;
  const busySource = Vue.ref(false);
  const busy = { status: 'busy', reason: 'Updating diagram. Retry after the update finishes.' };
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => { value = snapshot.value; };
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState,
    makeRef: Vue.ref,
    stepAvailability: (step) => (step.type === 'intent' && busySource.value ? busy : null)
  });
  const canUndoHistory = Vue.computed(() => { void history.revision.value; return history.canUndo(); });
  const canRedoHistory = Vue.computed(() => { void history.revision.value; return history.canRedo(); });
  await history.captureBaseline('initial');
  await history.runUndoable('Increment', () => { value = 1; });
  assert.equal(canUndoHistory.value, true);
  busySource.value = true;
  assert.equal(history.canUndo(), false);
  assert.equal(history.undoLabel(), '');
  assert.equal(canUndoHistory.value, false);
  assert.deepEqual(await history.undo(), busy);
  assert.equal(history.getUndoCount(), 1);
  busySource.value = false;
  assert.equal(canUndoHistory.value, true);
  assert.equal(history.undoLabel(), 'Increment');
  assert.equal(await history.undo(), true);
  assert.equal(canRedoHistory.value, true);
  busySource.value = true;
  assert.equal(history.canRedo(), false);
  assert.equal(history.redoLabel(), '');
  assert.equal(canRedoHistory.value, false);
  busySource.value = false;
  assert.equal(canRedoHistory.value, true);
  assert.equal(history.redoLabel(), 'Increment');
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
  const limitMessages = [];
  const history = createHistoryManager({
    maxBytes: 100,
    onHistoryLimit: (message) => limitMessages.push(message),
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
  // OV-1260: the Generate just taken stays undoable although it alone is over
  // the limit; it drops no older step, so nothing is announced.
  assert.equal(history.canUndo(), true);
  assert.deepEqual(limitMessages, []);
  assert.equal(await history.undo(), true);
  assert.equal(artifact.identity.compactSignature, 'a');
  assert.equal(history.canRedo(), true, 'the step just undone stays redoable');
  assert.equal(await history.redo(), true);
  assert.equal(artifact.identity.compactSignature, 'b');
}

{
  let artifact = {
    id: 'before',
    identity: { fingerprint: 'a'.repeat(64), compactSignature: 'before' },
    retainedBytes: 128,
    fileIds: []
  };
  let targetDraft = { startCoordinate: 1, reverseComplementOverride: false, anchorIntent: null };
  const history = createHistoryManager({
    buildIntent: async () => ({ unrelated: 'pending', targetDraft }),
    applyIntent: async (intent) => { targetDraft = intent.targetDraft; },
    buildCheckpoint: () => ({}),
    applyCheckpoint: async () => {},
    captureGeneratedArtifactHandle: () => artifact,
    restoreGeneratedArtifactHandle: async (handle) => { artifact = handle; },
    compareGeneratedArtifactHandles: (left, right) => (
      left.identity.fingerprint === right.identity.fingerprint
      && left.identity.compactSignature === right.identity.compactSignature
    )
  });
  const checkpointOptions = {
    shouldCommit: (result) => result.status === 'ok',
    captureIntentCheckpoint: () => structuredClone(targetDraft),
    restoreIntentCheckpoint: (checkpoint) => { targetDraft = structuredClone(checkpoint); }
  };
  await history.initializeIntentBaseline('Atomic artifact baseline');
  await history.runUndoableArtifactReplacement('Rotate record to feature', async () => {
    artifact = {
      id: 'after',
      identity: { fingerprint: 'b'.repeat(64), compactSignature: 'after' },
      retainedBytes: 128,
      fileIds: []
    };
    targetDraft = {
      startCoordinate: 75,
      reverseComplementOverride: true,
      anchorIntent: { schema: 1, recordKey: 'record-2' }
    };
    return { status: 'ok' };
  }, checkpointOptions);
  assert.equal(history.getUndoCount(), 1);
  await history.undo();
  assert.equal(artifact.id, 'before');
  assert.deepEqual(targetDraft, {
    startCoordinate: 1, reverseComplementOverride: false, anchorIntent: null
  });
  await history.redo();
  assert.equal(artifact.id, 'after');
  assert.equal(targetDraft.startCoordinate, 75);
  assert.equal(targetDraft.reverseComplementOverride, true);

  const undoCount = history.getUndoCount();
  const failed = await history.runUndoableArtifactReplacement(
    'Rejected record rotation',
    async (before) => {
      artifact = before;
      return { status: 'stale' };
    },
    checkpointOptions
  );
  assert.deepEqual(failed, { status: 'stale' });
  assert.equal(targetDraft.startCoordinate, 75);
  assert.equal(history.getUndoCount(), undoCount);
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
  const limitMessages = [];
  const history = createHistoryManager({
    maxActions: 2,
    onHistoryLimit: (message) => limitMessages.push(message),
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
  // OV-1365: the step-count window drops its oldest step silently; only the
  // byte limit is announced.
  assert.deepEqual(limitMessages, []);
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
  const limitMessages = [];
  const history = createHistoryManager({
    maxBytes: 220,
    onHistoryLimit: (message) => limitMessages.push(message),
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  await history.runUndoable('Small edit', () => {
    value = 'tiny';
  });
  assert.equal(history.getUndoCount(), 1);
  await history.runUndoable('Large state', () => {
    value = 'x'.repeat(1000);
  });
  // OV-1260: the newest step stays; the older one drops and the drop is announced.
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.undoLabel(), 'Large state');
  assert.deepEqual(limitMessages, ['Older undo history was discarded to stay within the history limit.']);
  assert.equal(await history.undo(), true);
  assert.equal(value, 'tiny');
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
    featureOverrides: {},
    featureStrokeOverrides: {},
    labelTextBulkOverrides: {},
    orthogroups: ref([]),
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    lastRunInfo: ref(null),
    pairwiseMatchFactors: ref({}),
    skipCaptureBaseConfig: ref(false),
    skipExtractOnSvgChange: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state: withTwoDrawings(state),
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
    applyConfigData: (_drawing, config) => {
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

  // Each drawing's checkpoint holds its own slots' files: Circular's inputs, and
  // Linear's comparison files and generated comparisons.
  const snapshot = await snapshots.drawingHistory('circular').buildCheckpoint();
  state.mode.value = 'linear';
  const linearSnapshot = await snapshots.drawingHistory('linear').buildCheckpoint();
  state.mode.value = 'circular';
  assert.equal(Object.hasOwn(snapshot.files, 'linearComparisons'), false);
  assert.equal(Object.hasOwn(linearSnapshot.files, 'c_gb'), false);
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
  await snapshots.drawingHistory('circular').applyCheckpoint(snapshot);
  state.mode.value = 'linear';
  await snapshots.drawingHistory('linear').applyCheckpoint(linearSnapshot);
  state.mode.value = 'circular';
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

  // The comparison files are the Linear drawing's: its History records them.
  state.mode.value = 'linear';
  const history = createHistoryManager({
    fileStore,
    ...snapshots.drawingHistory('linear')
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
  const serializedAfterConfigUndo = await snapshots.drawingHistory('linear').buildCheckpoint();
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
    featureOverrides: {},
    featureStrokeOverrides: {},
    labelTextBulkOverrides: {},
    orthogroups: ref([]),
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    lastRunInfo: ref(null),
    pairwiseMatchFactors: ref({}),
    skipCaptureBaseConfig: ref(false),
    skipExtractOnSvgChange: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state: withDrawings(state),
    fileStore,
    nextTick: async () => {
      if (!pendingModeReset) return;
      pendingModeReset = false;
      state.extractedFeatures.value = [];
      state.featureRecordIds.value = [];
      state.orthogroups.value = [];
      state.selectedOrthogroupId.value = '';
    },
    applyConfigData: (_drawing, config) => {
      state.form = { ...state.form, ...config.form };
      state.adv = { ...state.adv, ...config.adv };
    }
  });
  // PD-OI-090: a checkpoint restores its own drawing, the shown one; the
  // `mode` it holds is not written (a drawing switch is no History step).
  await snapshots.drawingHistory('circular').applyCheckpoint({
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

  assert.equal(state.mode.value, 'circular');
  assert.equal(pendingModeReset, false);
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
    featureOverrides: {},
    featureStrokeOverrides: {},
    labelTextBulkOverrides: {},
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
  state.legendEntries = Object.create({ value: [{ caption: 'tRNA', color: '#e8b441' }] });
  state.deletedLegendEntries = Object.create({ value: [] });
  const snapshots = createHistorySnapshotService({
    state: withDrawings(state),
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
  // R13: the composition owner is created after the snapshot service and
  // registers its capture once it exists. Until then the intent holds no
  // capture, as a service without that owner. U3b: the Legend list is intent,
  // so History captures no Legend row owners.
  const unregisteredIntent = await snapshots.drawingHistory('circular').buildIntent();
  assert.equal(Object.hasOwn(unregisteredIntent.ui, 'compositionUserDeltas'), false);
  assert.throws(() => snapshots.registerCapture('legendOrder', () => []), /Unknown History intent capture: legendOrder/);
  assert.throws(() => snapshots.registerCapture('legend', () => []), /Unknown History intent capture: legend/);
  const compositionRecord = { 'result-2': { primary: [[4, 5]] } };
  snapshots.registerCapture('composition', () => compositionRecord);
  const intent = await snapshots.drawingHistory('circular').buildIntent();
  assert.deepEqual(intent.editorState.legend.entries, [{ caption: 'tRNA', color: '#e8b441' }]);
  assert.deepEqual(intent.ui.compositionUserDeltas, compositionRecord);
  compositionRecord['result-2'].primary[0][0] = 9;
  assert.deepEqual(intent.ui.compositionUserDeltas['result-2'].primary, [[4, 5]], 'the intent holds a copy of the capture');
  snapshots.registerCapture('composition', null);
  assert.equal(Object.hasOwn((await snapshots.drawingHistory('circular').buildIntent()).ui, 'compositionUserDeltas'), false);
  state.legendEntries.value = [{ caption: 'tRNA', color: '#c026d3' }];
  await snapshots.drawingHistory('circular').applyIntent(intent, { changes: [{ path: ['editorState'] }] });
  assert.deepEqual(state.legendEntries.value, [{ caption: 'tRNA', color: '#e8b441' }]);
  assert.equal(Object.hasOwn(intent.editorState.legend, 'entryOwners'), false);
  assert.equal(forbiddenArtifactBuilds, 0);
  assert.equal(Object.prototype.hasOwnProperty.call(intent, 'results'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent, 'runState'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.features, 'extractedFeatures'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.editorState, 'featureCatalog'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.orthogroupState, 'groups'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.files, 'linearCanonicalComparisons'), false);
  // B23: a LOSAT-cache replay's comparison sequences are held by reference only.
  assert.deepEqual(intent.files.c_conservation_blasts, []);
  assert.deepEqual(
    Object.keys(intent.files.c_conservation_sequence_sources[0]).sort(),
    ['fileId', 'lastModified', 'name', 'size', 'type']
  );

  // B21: with uploaded BLAST rows, the optional comparison sequence is a user
  // choice, so the intent holds it and Undo restores it.
  state.files.c_conservation_blasts_source = null;
  const uploadIntent = await snapshots.drawingHistory('circular').buildIntent();
  assert.equal(uploadIntent.files.c_conservation_sequence_sources.length, 1);
  state.files.c_conservation_sequence_sources = [];
  await snapshots.drawingHistory('circular').applyIntent(uploadIntent, { changes: [{ path: ['files'] }] });
  assert.deepEqual(state.files.c_conservation_sequence_sources.map((file) => file?.name), ['generated-source.fa']);

  // B23: in a LOSAT-cache replay the ring rows and managed track slots name the
  // derived BLAST rows, so Undo and Redo of a ring removal restore those rows
  // and their comparison sequences with the ring files.
  const derivedBlasts = ['b', 'c', 'd'].map((id) => makeFile(`comparison-${id}.losatn.tsv`, 60));
  const derivedSources = ['b', 'c', 'd'].map((id) => makeFile(`comparison-${id}.fasta`, 1234));
  const ringFastas = ['b', 'c', 'd'].map((id) => makeFile(`comparison-${id}.fasta`, 1234));
  Object.assign(state.files, {
    c_conservation_blasts: [...derivedBlasts],
    c_conservation_blasts_source: 'losat-cache',
    c_conservation_fastas: [...ringFastas],
    c_conservation_sequence_sources: [...derivedSources]
  });
  const replayHistory = createHistoryManager({
    fileStore,
    collectCurrentFileIds: snapshots.collectCurrentFileIds,
    ...snapshots.drawingHistory('circular'),
    buildCheckpoint: () => {
      throw new Error('a ring removal must not build an artifact checkpoint');
    },
    applyCheckpoint: snapshots.drawingHistory('circular').applyCheckpoint
  });
  await replayHistory.initializeIntentBaseline('Loaded session');
  const replayRows = () => ({
    blasts: [...state.files.c_conservation_blasts],
    source: state.files.c_conservation_blasts_source,
    fastas: [...state.files.c_conservation_fastas],
    sequenceSources: [...state.files.c_conservation_sequence_sources]
  });
  const loadedRows = replayRows();
  await replayHistory.runUndoable('Remove series', () => {
    state.files.c_conservation_blasts = derivedBlasts.slice(1);
    state.files.c_conservation_fastas = ringFastas.slice(1);
    state.files.c_conservation_sequence_sources = derivedSources.slice(1);
  });
  const removedRows = replayRows();
  await replayHistory.undo();
  assert.deepEqual(replayRows(), loadedRows);
  replayRows().blasts.forEach((file, index) => assert.equal(file, derivedBlasts[index]));
  await replayHistory.redo();
  assert.deepEqual(replayRows(), removedRows);
  await replayHistory.undo();
  assert.deepEqual(replayRows(), loadedRows);
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
  const snapshots = createHistorySnapshotService({ state: withDrawings(state), fileStore });
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
    applyCheckpoint: snapshots.drawingHistory('circular').applyCheckpoint
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
    semanticFileWatchersSuppressed: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state: withDrawings(state),
    fileStore,
    buildUiStateData: () => ({}),
    applyUiStateData: () => {},
    serializeResults: () => [],
    applyResultsData: () => {}
  });
  // The drafts are Session-level values: the router's project stream records them (R15-23).
  const history = createHistoryRouter({
    drawingIds: ['circular'],
    activeDrawingId: () => 'circular',
    managerOptions: () => ({ fileStore, ...snapshots.drawingHistory('circular') }),
    sessionHistory: snapshots.sessionHistory()
  });
  const draftState = () => ({
    selectedAnnotation: structuredClone(state.selectedAnnotation.value),
    selectedSpecificPreset: state.selectedSpecificPreset.value,
    newSpecRule: structuredClone(state.newSpecRule),
    newPriorityRule: structuredClone(state.newPriorityRule),
    newColorFeat: state.newColorFeat.value,
    newColorVal: state.newColorVal.value,
    newFeatureToAdd: state.newFeatureToAdd.value
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
  });
  await history.undo();
  assert.deepEqual(draftState(), after);

  // R15-2 retired Add legend item and its form drafts. A Session-level step
  // that still carries them restores the other drafts and writes nothing for
  // them, even where a holder of the old name exists. A drawing's checkpoint
  // holds no drafts and writes none (R15-23).
  const retiredDrafts = { newLegendCaption: 'Old draft', newLegendColor: '#123123' };
  const sessionPorts = snapshots.sessionHistory();
  const oldSession = { ...sessionPorts.buildIntent() };
  oldSession.drafts = { ...oldSession.drafts, ...retiredDrafts };
  const oldCheckpoint = { ...(await snapshots.drawingHistory('circular').buildCheckpoint()) };
  assert.equal(Object.hasOwn(oldCheckpoint, 'drafts'), false);
  oldCheckpoint.drafts = { ...oldSession.drafts, newColorVal: '#000003' };
  state.newLegendCaption = ref('');
  state.newLegendColor = ref('#808080');
  state.newColorVal.value = '#000001';
  await sessionPorts.applyIntent(oldSession);
  assert.deepEqual(draftState(), after);
  state.newColorVal.value = '#000002';
  await snapshots.drawingHistory('circular').applyCheckpoint(oldCheckpoint);
  assert.equal(state.newColorVal.value, '#000002');
  assert.deepEqual([state.newLegendCaption.value, state.newLegendColor.value], ['', '#808080']);
  assert.equal(Object.hasOwn(sessionPorts.buildIntent().drafts, 'newLegendCaption'), false);
  assert.equal(Object.hasOwn(await snapshots.drawingHistory('circular').buildIntent(), 'drafts'), false);
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
  const noticesA = [{
    recordKey: 'record-a', biologicalFeatureId: 'gone', status: 'unresolved', kinds: ['label_text'], resultIndex: 0
  }];
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
    featureRecordIds: ref(['record-a']),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: { 'feature-a': '#112233' },
    featureVisibilityManualRules: [],
    featureOverrides: {},
    // Python's identity notices are an artifact-owned value (design Q4 6.3).
    featureIdentityNotices: ref(noticesA),
    canonicalLabelOverrideRows: ref([]),
    labelTextBulkOverrides: {},
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
  // D-37: installing a Generate's artifact over other shown Results clears the
  // similarity-group copy feedback, also when both share one groups array.
  let copyFeedbackClears = 0;
  const snapshots = createHistorySnapshotService({
    state: withDrawings(state),
    fileStore,
    nextTick: async () => {},
    clearSimilarityGroupCopyFeedback: () => { copyFeedbackClears += 1; },
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
    applyUiStateData: (_drawing, ui) => {
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
    applyFeatureStateData: (_drawing, features) => {
      state.extractedFeatures.value = features.extractedFeatures;
      state.biologicalFeatures.value = features.biologicalFeatures;
      state.featureRecordIds.value = features.featureRecordIds;
      state.selectedFeatureRecordIdx.value = features.selectedFeatureRecordIdx;
      Object.assign(state.featureColorOverrides, features.featureColorOverrides);
    },
    applyEditorStateData: (_drawing, editorState, options) => {
      assert.equal(options.normalized, true);
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

  const handleA = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
  {
    // A reused Generate's commit: new Results with the shown groups array.
    const clears = copyFeedbackClears;
    const reusedGenerate = { ...handleA.ownerSet, results: [{ ...resultA }] };
    snapshots.installGeneratedArtifactOwnerSet(reusedGenerate);
    assert.equal(state.orthogroups.value, groupsA);
    assert.equal(copyFeedbackClears, clears + 1);
    // Installing the shown owner set again (a rollback before activation): no clear.
    snapshots.installGeneratedArtifactOwnerSet(reusedGenerate);
    assert.equal(copyFeedbackClears, clears + 1);
    snapshots.installGeneratedArtifactOwnerSet(handleA.ownerSet);
    assert.equal(copyFeedbackClears, clears + 2);
    assert.equal(state.results.value, handleA.ownerSet.results);
  }
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
  assert.equal(handleA.ownerSet.featureIdentityNotices, noticesA);
  assert.ok(handleA.retainedBytes >= 400_000);

  state.resultGenerationKey.value += 1;
  state.zoom.value = 2;
  state.canvasPan.x = 120;
  const transientOnlyHandle = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
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
  const poisonNotices = new Proxy([], {
    get() { throw new Error('Identity notices must not be traversed during capture.'); }
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
  state.featureIdentityNotices.value = poisonNotices;
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
  const poisonHandle = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
  assert.equal(poisonHandle.ownerSet.results[0], poisonResult);
  assert.equal(poisonHandle.ownerSet.featureCatalog, poisonCatalog);
  assert.equal(poisonHandle.ownerSet.losatCache, poisonCache);
  assert.equal(poisonHandle.ownerSet.matchSequenceOwner, poisonSequenceOwner);
  assert.equal(poisonHandle.ownerSet.featureIdentityNotices, poisonNotices);
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
  const invalidatedIdentityHandle = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
  assert.equal(
    invalidatedIdentityHandle.identity.fingerprint,
    '',
    'Replacing a Result object must invalidate its Worker-byte identity even at equal length.'
  );
  state.results.value = [{ name: 'edited.svg', content: '<svg id="edited" />' }];
  state.featureColorOverrides['feature-a'] = '#abcdef';
  state.featureIdentityNotices.value = [{ recordKey: 'record-b', biologicalFeatureId: 'bio-b' }];
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
  const handleB = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());

  assert.equal(handleA.ownerSet.results[0], resultA);
  assert.equal(handleA.mutableIntent.features.featureColorOverrides['feature-a'], '#112233');
  assert.equal(handleA.mutableIntent.editorState.legend.entries[0].caption, 'A');
  assert.equal(handleA.ownerSet.featureIdentityNotices, noticesA);
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
  const clearsBeforeUndo = copyFeedbackClears;
  await snapshots.restoreGeneratedArtifactHandle(handleA);
  // Undo of a Generate step clears the copy feedback.
  assert.equal(copyFeedbackClears, clearsBeforeUndo + 1);
  assert.equal(resultAdmissionCalls, 0);
  assert.equal(state.results.value[0], resultA);
  assert.equal(state.extractedFeatures.value, extractedA);
  assert.equal(state.biologicalFeatures.value, biologicalA);
  assert.equal(state.featureCatalog.value, catalogA);
  assert.equal(state.orthogroups.value, groupsA);
  assert.equal(state.featureColorOverrides['feature-a'], '#112233');
  assert.equal(state.legendEntries.value[0].caption, 'A');
  assert.equal(state.featureIdentityNotices.value, noticesA);
  assertAuthorityOwners(authorityOwnersA);
  assert.equal(state.appliedPaletteColors.value, appliedPaletteColorsA);
  assert.equal(state.pendingPaletteColors.value, pendingPaletteColorsA);
  const recapturedA = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
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
    captureGeneratedArtifactHandle: () => snapshots.captureGeneratedArtifactHandle(state.activeDrawing()),
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

  const restoredHandleBeforePaletteReplacement = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
  const equivalentAppliedPaletteColors = { ...state.appliedPaletteColors.value };
  const equivalentPendingPaletteColors = { ...state.pendingPaletteColors.value };
  state.appliedPaletteColors.value = equivalentAppliedPaletteColors;
  state.pendingPaletteColors.value = equivalentPendingPaletteColors;
  const equivalentPaletteHandle = snapshots.captureGeneratedArtifactHandle(state.activeDrawing());
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

{
  let file = 'original.tsv';
  let rules = ['original'];
  let svg = 'original legend';
  const history = createHistoryManager({
    buildIntent: () => ({ file, rules }),
    applyIntent: (intent) => { file = intent.file; rules = [...intent.rules]; },
    buildCheckpoint: () => ({ file, rules: [...rules], svg }),
    applyCheckpoint: (checkpoint) => {
      file = checkpoint.file;
      rules = [...checkpoint.rules];
      svg = checkpoint.svg;
    }
  });
  await history.initializeIntentBaseline('Rejected artifact baseline');
  await history.runUndoableCheckpoint('Import specific colors', async () => {
    file = 'colors.tsv';
    await history.runUndoableCheckpoint('Commit prepared rules', async () => {
      rules = ['canonical caption'];
      svg = 'canonical legend';
      await history.runUndoable('Derived preview', () => { svg += ' and colors'; });
    });
  });
  assert.equal(history.getUndoCount(), 1, 'an upload and its canonical rule/legend edit are one action');
  assert.equal(history.getDiagnostics().artifactCheckpointBuilds, 2, 'nested edits share the outer capture');
  await history.undo();
  assert.deepEqual({ file, rules, svg }, { file: 'original.tsv', rules: ['original'], svg: 'original legend' });
  await history.redo();
  assert.deepEqual({ file, rules, svg }, { file: 'colors.tsv', rules: ['canonical caption'], svg: 'canonical legend and colors' });
  await assert.rejects(history.runUndoableCheckpoint('Invalid import', async () => {
    throw new Error('invalid known value');
  }), /invalid known value/);
  assert.equal(history.getUndoCount(), 1);
  await history.runUndoable('Next edit', () => { rules = ['next']; });
  assert.equal(history.getUndoCount(), 2, 'a failed checkpoint releases its nesting scope');
  await history.undo();
  assert.deepEqual(rules, ['canonical caption']);
}

for (const previousSuppressed of [false, true]) {
  const state = { semanticFileWatchersSuppressed: ref(previousSuppressed) };
  const snapshots = createHistorySnapshotService({
    state: withDrawings(state),
    fileStore: createHistoryFileStore(),
    nextTick: async () => {
      assert.equal(state.semanticFileWatchersSuppressed.value, true);
      throw new Error('restore did not settle');
    }
  });
  await assert.rejects(snapshots.drawingHistory('circular').applyCheckpoint({}), /restore did not settle/);
  assert.equal(state.semanticFileWatchersSuppressed.value, previousSuppressed,
    'checkpoint failure restores the prior semantic file-watcher suppression');
}

{
  const before = { identity: { fingerprint: '', compactSignature: 'before' }, retainedBytes: 0, fileIds: [] };
  const primary = {code:'RESULT_INVALID',stage:'result-admission',message:'PRIVATE_PRIMARY_SENTINEL'};
  const history = createHistoryManager({buildIntent:async()=>({}),applyIntent:async()=>{},
    buildCheckpoint:async()=>({}),applyCheckpoint:async()=>{},
    captureGeneratedArtifactHandle:()=>before,restoreGeneratedArtifactHandle:async()=>{},
    compareGeneratedArtifactHandles:()=>false});
  await history.initializeIntentBaseline('Rejected artifact baseline');
  let rollbackSettled = false;
  await assert.rejects(history.runUndoableArtifactReplacement('Rejected artifact', async()=>{throw primary;}, {
    restoreAppliedArtifact:async()=>{await Promise.resolve();rollbackSettled=true;throw new Error('PRIVATE_ROLLBACK_SENTINEL');}
  }), error=>error===primary && error.artifactRestoreFailed===true);
  assert.equal(rollbackSettled,true);
  assert.equal(history.getUndoCount(),0);
  assert.equal(history.getRedoCount(),0);
}

// OV-460: an undoable edit made while an async undoable action is still
// running applies at once and is its own step under its own name. It, a drag
// end, or a settlement commits the action's transaction while its `fn` still
// awaits; what `fn` writes afterwards is recorded as a step under the action's
// label, so Undo then Redo keeps it.
{
  const START = { label: 'alpha', table: null, color: 'red' };
  const hold = () => {
    let release;
    const held = new Promise((resolve) => { release = resolve; });
    return { held, release };
  };
  const tick = () => new Promise((resolve) => setImmediate(resolve));
  /** @param {(harness: Record<string, any>) => Record<string, any>} [owners] More History owners (the artifact owner). */
  const createActionHarness = async (owners = () => ({})) => {
    const harness = { intent: structuredClone(START), restoreHeld: null };
    harness.history = createHistoryManager({
      buildIntent: () => structuredClone(harness.intent),
      applyIntent: async (restored) => {
        await harness.restoreHeld;
        harness.intent = structuredClone(restored);
      },
      buildCheckpoint: () => structuredClone(harness.intent),
      applyCheckpoint: (checkpoint) => { harness.intent = structuredClone(checkpoint); },
      ...owners(harness)
    });
    await harness.history.initializeIntentBaseline();
    harness.undoAll = async () => {
      const labels = [];
      while (harness.history.canUndo()) {
        labels.push(harness.history.undoLabel());
        await harness.history.undo();
      }
      return labels;
    };
    harness.redoAll = async () => { while (harness.history.canRedo()) await harness.history.redo(); };
    // Load Label TSV while Python evaluates the table: the action awaits, then writes.
    harness.heldLoad = (write = () => { harness.intent.table = 'applied'; }) => {
      const { held, release } = hold();
      const done = harness.history.runUndoable('Load label edits', async () => {
        await held;
        write();
      });
      return { done, release };
    };
    return harness;
  };

  // A second edit joins the open transaction and commits it under its own
  // name; the load then writes.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    load.release();
    await load.done;
    assert.equal(h.history.mutationPending(), false);
    assert.equal(h.history.getUndoCount(), 2, 'the later write is a step');
    assert.deepEqual(await h.undoAll(), ['Load label edits', 'Change label text'],
      'the edit and the action\'s later write are each a step under their own label');
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'EDITED', table: 'applied', color: 'red' }, 'Redo keeps both writes');
  }

  // The running action ends first: its step keeps its name, and what the
  // later edit writes afterwards is a step under the edit's name.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    const edit = hold();
    await tick();
    const recolor = h.history.runUndoable('Change feature color', async () => {
      await edit.held;
      h.intent.color = 'blue';
    });
    load.release();
    await load.done;
    edit.release();
    await recolor;
    assert.deepEqual(await h.undoAll(), ['Change feature color', 'Load label edits']);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'alpha', table: 'applied', color: 'blue' });
  }

  // A call the running action makes itself keeps the step and its name: one
  // made synchronously inside `fn`, and one made after an await that says so
  // (`joinsOpenStep`: the rule table commit of a color action, a dialog choice).
  {
    const h = await createActionHarness();
    await h.history.runUndoable('Change feature color', () => {
      h.intent.color = 'blue';
      return h.history.runUndoable('Change specific color rules', () => { h.intent.table = 'rule'; });
    });
    assert.deepEqual(await h.undoAll(), ['Change feature color']);
    assert.deepEqual(h.intent, START);

    const prepared = hold();
    const recolor = h.history.runUndoable('Change feature color', async () => {
      h.intent.color = 'blue';
      await prepared.held;
      await h.history.runUndoable('Change specific color rules', () => { h.intent.table = 'rule'; }, { joinsOpenStep: true });
    });
    await tick();
    prepared.release();
    await recolor;
    assert.deepEqual(await h.undoAll(), ['Change feature color']);
    assert.deepEqual(h.intent, START);
  }

  // The same color action made while the load runs: its own call and the
  // call that joins it commit one step under the color action's name.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    const prepared = hold();
    await tick();
    const recolor = h.history.runUndoable('Change feature color', async () => {
      h.intent.color = 'blue';
      await prepared.held;
      await h.history.runUndoable('Change specific color rules', () => { h.intent.label = 'rule'; }, { joinsOpenStep: true });
    });
    await tick();
    prepared.release();
    await recolor;
    load.release();
    await load.done;
    assert.deepEqual(await h.undoAll(), ['Load label edits', 'Change feature color']);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'rule', table: 'applied', color: 'blue' });
  }

  // R11: an action that joins its control's transaction keeps the control's
  // name when no action runs on that transaction, also while an action whose
  // step was already committed is still running.
  {
    const h = await createActionHarness();
    const button = async () => h.history.begin('Change setting', { source: 'input-adapter', owner: { tagName: 'BUTTON' } });
    await button();
    await h.history.runUndoable('Change feature color', () => { h.intent.color = 'blue'; });
    assert.deepEqual(await h.undoAll(), ['Change setting']);

    const load = h.heldLoad();
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    await button();
    await h.history.runUndoable('Change feature color', () => { h.intent.color = 'blue'; });
    load.release();
    await load.done;
    assert.deepEqual(await h.undoAll(), ['Load label edits', 'Change setting', 'Change label text']);
    assert.deepEqual(h.intent, START);
  }

  // The same when the load writes nothing afterwards: no step and no bookkeeping.
  {
    const h = await createActionHarness();
    const load = h.heldLoad(() => {});
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    const revision = h.history.revision.value;
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 1);
    assert.equal(h.history.revision.value, revision);
    await h.undoAll();
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'EDITED', table: null, color: 'red' });
  }

  // OV-461: a drag joins the running action's transaction and commits it at
  // the drag end; the action's later write is no longer lost.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    await tick();
    const drag = await h.history.begin('Move legend', { source: 'legend-drag', owner: Symbol('Move legend') });
    h.intent.color = 'blue';
    await h.history.commit(drag);
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 2);
    assert.equal(h.history.undoLabel(), 'Load label edits');
    await h.undoAll();
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'alpha', table: 'applied', color: 'blue' });
  }

  // The remainder is recorded while a control's transaction is open: that
  // transaction continues from the remainder, so no change is in two steps.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    const typing = await h.history.begin('Change setting', { source: 'input-adapter', owner: { tagName: 'INPUT' } });
    h.intent.color = 'bl';
    load.release();
    await load.done;
    h.intent.color = 'blue';
    await h.history.commit(typing);
    assert.equal(h.history.getUndoCount(), 3);
    assert.equal(h.history.undoLabel(), 'Change setting');
    await h.history.undo();
    assert.deepEqual(h.intent, { label: 'EDITED', table: 'applied', color: 'bl' }, 'the control step holds only what followed the remainder');
    assert.equal(h.history.undoLabel(), 'Load label edits');
    await h.history.undo();
    assert.deepEqual(h.intent, { label: 'EDITED', table: null, color: 'red' });
    await h.undoAll();
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'EDITED', table: 'applied', color: 'blue' });
  }

  // A cancelled transaction has no remainder (a Session replacement resets the baseline).
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    await tick();
    await h.history.initializeIntentBaseline();
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0);
    assert.equal(h.history.mutationPending(), false);
  }

  // A checkpoint settles the load's transaction and is still open when the
  // load writes: the checkpoint's step holds the write, as before.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    const imported = hold();
    await tick();
    const colors = h.history.runUndoableCheckpoint('Import specific colors', async () => {
      h.intent.color = 'blue';
      await imported.held;
    });
    await tick();
    load.release();
    await load.done;
    imported.release();
    await colors;
    assert.deepEqual(await h.undoAll(), ['Import specific colors']);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'alpha', table: 'applied', color: 'blue' });
  }

  // An Undo that is being applied owns the stacks: the load ends meanwhile
  // and records nothing, so the Undo pops the step it restored.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    const restore = hold();
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    h.restoreHeld = restore.held;
    const undo = h.history.undo();
    await tick();
    load.release();
    await load.done;
    restore.release();
    await undo;
    h.restoreHeld = null;
    assert.deepEqual([h.history.getUndoCount(), h.history.getRedoCount()], [0, 1]);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'EDITED', table: null, color: 'red' });
  }

  // OV-463: the action's later write is recorded from the last recorded
  // intent, which a Generate does not refresh and a Session Load replaces.
  // After either, or while a Generate is open, the write is in no step.
  const artifactOf = (id) => ({ id, identity: { fingerprint: id, compactSignature: id }, retainedBytes: 8, fileIds: [] });
  // Generate replaces the Result, and the Legend entries follow the Result.
  const withResult = (harness) => {
    harness.artifact = artifactOf('a');
    harness.generate = (held = null) => harness.history.runUndoableArtifactReplacement('Generate diagram', async () => {
      harness.artifact = artifactOf('b');
      harness.intent.color = 'legend-of-b';
      await held;
      return { status: 'ok' };
    }, { shouldCommit: (result) => result.status === 'ok' });
    return {
      captureGeneratedArtifactHandle: () => harness.artifact,
      restoreGeneratedArtifactHandle: async (handle) => {
        harness.artifact = handle;
        harness.intent.color = handle.id === 'a' ? 'red' : 'legend-of-b';
      }
    };
  };
  for (const write of [() => {}, undefined]) {
    // Generate settles the load's step and ends; the load then ends, declined
    // (it writes nothing) or not.
    const h = await createActionHarness(withResult);
    const load = h.heldLoad(write);
    await tick();
    await h.generate();
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 1, 'Generate\'s changes are not a step of the load');
    assert.equal(h.history.undoLabel(), 'Generate diagram');
    await h.history.undo();
    assert.deepEqual([h.artifact.id, h.intent.color], ['a', 'red'], 'one Undo returns the Result and what follows it');
    await h.history.redo();
    assert.deepEqual([h.artifact.id, h.intent.color], ['b', 'legend-of-b']);

    // An action that starts after the Generate records its later write again.
    const later = h.heldLoad(() => { h.intent.table = 'later'; });
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    later.release();
    await later.done;
    assert.equal(h.history.getUndoCount(), 3);
    assert.equal(h.history.undoLabel(), 'Load label edits');
    await h.history.undo();
    assert.deepEqual(h.intent, { label: 'EDITED', table: write ? null : 'applied', color: 'legend-of-b' },
      'the step holds only what that action wrote');
  }
  {
    // The load ends while the Generate is still open.
    const h = await createActionHarness(withResult);
    const load = h.heldLoad();
    const generating = hold();
    await tick();
    const generated = h.generate(generating.held);
    await tick();
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0, 'no step under the open Generate');
    generating.release();
    await generated;
    assert.deepEqual(await h.undoAll(), ['Generate diagram']);
  }
  for (const baseline of ['initializeIntentBaseline', 'captureBaseline']) {
    // A Session Load replaces the baseline after the load's step was
    // committed; what the stale load writes is not a step of the new Session.
    const h = await createActionHarness();
    const load = h.heldLoad();
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    h.intent = { label: 'session', table: null, color: 'session' };
    await h.history[baseline]('Loaded session');
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0);
    assert.equal(await h.history.hasChangesSinceSavePoint(), false);
  }

  // A checkpoint or an artifact replacement that ends without committing
  // leaves its writes in no step, and they are not the running action's.
  {
    const h = await createActionHarness();
    const load = h.heldLoad(() => {});
    await tick();
    await h.history.runUndoableCheckpoint('Change uploaded file', async () => { h.intent.table = 'declined'; return false; },
      { shouldCommit: (result) => result !== false });
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0, 'a declined checkpoint is in no step');
    assert.equal(await h.history.hasChangesSinceSavePoint(), false);
  }
  {
    const h = await createActionHarness();
    const load = h.heldLoad(() => {});
    await tick();
    await assert.rejects(h.history.runUndoableCheckpoint('Change uploaded file', async () => {
      h.intent.table = 'failed';
      throw new Error('import failed');
    }), /import failed/);
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0, 'a checkpoint that throws is in no step');
  }
  {
    const h = await createActionHarness(withResult);
    const load = h.heldLoad(() => {});
    await tick();
    await h.history.runUndoableArtifactReplacement('Generate diagram', async () => { h.intent.color = 'partial'; return { status: 'cancelled' }; },
      { shouldCommit: (result) => result.status === 'ok' });
    load.release();
    await load.done;
    assert.equal(h.history.getUndoCount(), 0, 'a Generate that is not committed is in no step');
  }

  // OV-464: a call that joins commits under the step name of the action it
  // belongs to (the one started on the step with its label), also when a
  // later action has renamed the step. Here: a popup recolor from a color
  // input, then Reset color while the recolor's rule write is prepared.
  {
    const h = await createActionHarness();
    const ruleAction = (label, write) => {
      const prepared = hold();
      const done = h.history.runUndoable(label, async () => {
        await prepared.held;
        await h.history.runUndoable(label, write, { joinsOpenStep: true });
      });
      return { done, release: prepared.release };
    };
    await h.history.begin('Change color', { source: 'input-adapter', owner: { tagName: 'INPUT' } });
    const recolor = ruleAction('Change feature color', () => { h.intent.color = 'blue'; });
    await tick();
    await h.history.begin('Change setting', { source: 'input-adapter', owner: { tagName: 'BUTTON' } });
    const reset = ruleAction('Reset feature color', () => { h.intent.color = 'red'; h.intent.table = 'reset'; });
    await tick();
    recolor.release();
    await recolor.done;
    reset.release();
    await reset.done;
    assert.equal(h.history.mutationPending(), false);
    assert.deepEqual(await h.undoAll(), ['Reset feature color', 'Change color'],
      'the recolor keeps the name its step started with');
    assert.deepEqual(h.intent, START);

    // The same two actions while a load runs: each is a step under its own
    // name, and so is what the load writes afterwards.
    const load = h.heldLoad();
    await tick();
    const first = ruleAction('Change feature color', () => { h.intent.color = 'blue'; });
    await tick();
    const second = ruleAction('Reset feature color', () => { h.intent.label = 'reset'; });
    await tick();
    first.release();
    await first.done;
    second.release();
    await second.done;
    load.release();
    await load.done;
    assert.deepEqual(await h.undoAll(), ['Load label edits', 'Reset feature color', 'Change feature color']);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'reset', table: 'applied', color: 'blue' });
  }

  // Two actions started with the same label on one step: a joining call
  // belongs to the first of them, whichever rule write arrives first.
  for (const order of [[0, 1], [1, 0]]) {
    const g = await createActionHarness();
    const sameLabelAction = (label, write) => {
      const prepared = hold();
      const done = g.history.runUndoable(label, async () => {
        await prepared.held;
        await g.history.runUndoable(label, write, { joinsOpenStep: true });
      });
      return { done, release: prepared.release };
    };
    await g.history.begin('Change color', { source: 'input-adapter', owner: { tagName: 'INPUT' } });
    const one = sameLabelAction('Change feature color', () => { g.intent.color = 'blue'; });
    await tick();
    await g.history.begin('Change color', { source: 'input-adapter', owner: { tagName: 'INPUT' } });
    const two = sameLabelAction('Change feature color', () => { g.intent.table = 'two'; });
    await tick();
    for (const index of order) {
      [one, two][index].release();
      await [one, two][index].done;
    }
    assert.equal(g.history.mutationPending(), false);
    assert.deepEqual(await g.undoAll(), ['Change feature color', 'Change color']);
  }

  // Three overlapping actions: every write is in a step under the name of
  // the action that committed it, and Undo then Redo keeps all of them.
  {
    const h = await createActionHarness();
    const load = h.heldLoad();
    const recoloring = hold();
    await tick();
    const recolor = h.history.runUndoable('Change feature color', async () => {
      h.intent.color = 'pre';
      await recoloring.held;
      h.intent.color = 'blue';
    });
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    recoloring.release();
    await recolor;
    load.release();
    await load.done;
    assert.deepEqual(await h.undoAll(), ['Load label edits', 'Change feature color', 'Change label text']);
    assert.deepEqual(h.intent, START);
    await h.redoAll();
    assert.deepEqual(h.intent, { label: 'EDITED', table: 'applied', color: 'blue' });
  }

  // An action that throws after another edit committed its step: what it
  // wrote before the throw is in no step, and no step stays open.
  {
    const h = await createActionHarness();
    const failing = hold();
    const load = h.history.runUndoable('Load label edits', async () => {
      await failing.held;
      h.intent.table = 'partial';
      throw new Error('load failed');
    });
    await tick();
    await h.history.runUndoable('Change label text', () => { h.intent.label = 'EDITED'; });
    failing.release();
    await assert.rejects(load, /load failed/);
    assert.equal(h.history.mutationPending(), false);
    assert.deepEqual(await h.undoAll(), ['Change label text']);
  }

  // Random interleavings of async actions (with and without a rule write
  // that joins), synchronous edits, typing, drags, commands, checkpoints,
  // Undo and Redo, from a fixed seed. Every action writes when it starts and
  // again when it ends. After all of them: no step is open, Undo of
  // everything gives the start, and Redo of those steps gives the end.
  {
    let seed = 20261012;
    const random = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
    const pick = (list) => list[Math.floor(random() * list.length)];
    const KEYS = ['a', 'b', 'c', 'd'];
    const FUZZ_START = { a: 0, b: 0, c: 0, d: 0, command: 0 };
    for (let run = 0; run < 600; run += 1) {
      const traverses = run % 2 === 1;
      let intent = structuredClone(FUZZ_START);
      let counter = 0;
      const log = [];
      const history = createHistoryManager({
        buildIntent: () => structuredClone(intent),
        applyIntent: async (restored) => { intent = structuredClone(restored); },
        buildCheckpoint: () => structuredClone(intent),
        applyCheckpoint: (checkpoint) => { intent = structuredClone(checkpoint); }
      });
      await history.initializeIntentBaseline();
      const write = (who) => { const key = pick(KEYS); intent[key] = ++counter; log.push(`${who} writes ${key}`); };
      const actions = [];
      const startAction = (joins) => {
        const id = `action ${actions.length + 1}`;
        const held = hold();
        log.push(`${id} starts`);
        const done = history.runUndoable(id, async () => {
          write(id);
          await held.held;
          if (joins) await history.runUndoable(id, () => write(`${id} rule`), { joinsOpenStep: true });
          else write(id);
        });
        actions.push({ id, done, release: held.release, ended: false });
      };
      const end = async (action) => { action.ended = true; log.push(`${action.id} ends`); action.release(); await action.done; };
      const controls = [{ id: 'input 1', tx: null }, { id: 'input 2', tx: null }];
      const type = async (control) => {
        if (!control.tx || control.tx.closed) control.tx = await history.begin(`Edit ${control.id}`, { source: 'input-adapter', owner: control });
        write(control.id);
      };
      // The input adapter commits on blur unless an action owns the commit.
      const blur = async (control) => {
        const tx = control.tx;
        if (!tx || tx.deferAdapterCommit) return;
        control.tx = null;
        if (!tx.closed) { log.push(`${control.id} commits`); await history.commit(tx); }
      };
      const kinds = ['start', 'startJoining', 'edit', 'type', 'type', 'blur', 'drag', 'command', 'checkpoint'];
      for (let step = 4 + Math.floor(random() * 8); step > 0; step -= 1) {
        const running = actions.filter((action) => !action.ended);
        const kind = pick([...kinds, ...(running.length ? ['end', 'end', 'end'] : []),
          ...(traverses ? ['undo', 'redo', 'heldCheckpoint'] : [])]);
        if ((kind === 'start' || kind === 'startJoining') && actions.length < 3) startAction(kind === 'startJoining');
        else if (kind === 'edit') { log.push('edit'); await history.runUndoable('edit', () => write('edit')); }
        else if (kind === 'type') await type(pick(controls));
        else if (kind === 'blur') await blur(pick(controls));
        else if (kind === 'drag') {
          const tx = await history.begin('drag', { owner: Symbol('drag') });
          write('drag');
          await history.commit(tx);
        } else if (kind === 'command') {
          const before = intent.command;
          const after = ++counter;
          log.push('command');
          await history.runUndoableCommand('command', () => ({
            apply: () => { intent.command = after; }, revert: () => { intent.command = before; }
          }));
        } else if (kind === 'checkpoint') await history.runUndoableCheckpoint('checkpoint', () => write('checkpoint'));
        else if (kind === 'undo' && history.canUndo()) { log.push('undo'); await history.undo(); }
        else if (kind === 'redo' && history.canRedo()) { log.push('redo'); await history.redo(); }
        else if (kind === 'heldCheckpoint' && running.length) {
          // An action ends inside a checkpoint that is still open.
          const held = hold();
          const checkpoint = history.runUndoableCheckpoint('held checkpoint', async () => { write('held checkpoint'); await held.held; });
          await tick();
          await end(pick(running));
          held.release();
          await checkpoint;
        } else if (kind === 'end' && running.length) await end(pick(running));
        await tick();
      }
      for (const action of actions) if (!action.ended) await end(action);
      for (const control of controls) await blur(control);
      await tick();
      const trace = `run ${run}: ${log.join('; ')}`;
      // A control transaction an ended action deferred is settled by the next owner.
      if (history.mutationPending()) await history.commit(await history.begin('settle', { owner: Symbol('settle') }));
      assert.equal(history.mutationPending() || Boolean(history.historyAvailability()), false, trace);
      const ended = structuredClone(intent);
      let undone = 0;
      while (history.canUndo()) { await history.undo(); undone += 1; }
      assert.deepEqual(intent, FUZZ_START, trace);
      for (; undone > 0; undone -= 1) await history.redo();
      assert.deepEqual(intent, ended, trace);
    }
  }
  console.log('history action remainder tests passed');
}

// OV-1260: a kept rerender that grows the last step over the byte limit
// drops the Redo steps, and the drop is announced as an Undo drop is.
{
  let rules = 0;
  let figure = 'drawn-0';
  const limitMessages = [];
  const history = createHistoryManager({
    maxBytes: 1000,
    onHistoryLimit: (message) => limitMessages.push(message),
    buildIntent: () => ({ rules }),
    applyIntent: (intent) => { rules = intent.rules; },
    buildCheckpoint: () => ({ rules, figure }),
    applyCheckpoint: (checkpoint) => { rules = checkpoint.rules; figure = checkpoint.figure; }
  });
  await history.captureBaseline('initial');
  await history.runUndoableCheckpoint('Change feature color', () => {
    rules = 1;
    figure = 'shown-1';
  }, { drawn: { before: () => true, change: false } });
  await history.runUndoable('Edit setting', () => { rules = 2; });
  await history.undo();
  assert.equal(history.getRedoCount(), 1);
  assert.deepEqual(limitMessages, []);
  figure = 'x'.repeat(2000);
  assert.equal(await history.keepRedrawnResult(), true);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(history.undoLabel(), 'Change feature color');
  assert.deepEqual(limitMessages, ['Redo history was discarded to stay within the history limit.']);
}

// D-24 (OV-731): a checkpoint step's caller says whether the figure of each
// side is the one its state draws, and Undo and Redo pass that on. A side not
// yet drawn takes the figure of a rerender drawn from the current state (D-36)
// only while its step is the last one and nothing is open.
{
  let rules = 0;
  let figure = 'drawn-0';
  const restores = [];
  const history = createHistoryManager({
    buildIntent: () => ({ rules }),
    applyIntent: (intent) => { rules = intent.rules; },
    buildCheckpoint: () => ({ rules, figure }),
    applyCheckpoint: (checkpoint, restored) => {
      rules = checkpoint.rules;
      figure = checkpoint.figure;
      restores.push(restored);
    }
  });
  await history.captureBaseline('initial');
  await history.runUndoableCheckpoint('Change feature color', () => {
    rules = 1;
    figure = 'shown-1';
  }, { drawn: { before: () => true, change: false } });
  figure = 'drawn-1';
  assert.equal(await history.keepRedrawnResult(), true);
  assert.equal(await history.keepRedrawnResult(), false, 'the kept side is drawn');
  await history.undo();
  await history.redo();
  assert.deepEqual({ rules, figure }, { rules: 1, figure: 'drawn-1' });
  assert.deepEqual(restores, [{ drawn: true }, { drawn: true }]);

  // A later step, or an Undo of the step, keeps the figure out of it.
  await history.runUndoableCheckpoint('Change feature color again', () => {
    rules = 2;
    figure = 'shown-2';
  }, { drawn: { before: () => false, change: true } });
  await history.runUndoable('Edit setting', () => { rules = 3; });
  figure = 'drawn-3';
  assert.equal(await history.keepRedrawnResult(), false, 'not a checkpoint step');
  await history.undo();
  await history.undo();
  assert.deepEqual(restores.at(-1), { drawn: false }, 'the before side was not drawn');
  // The step before is the last one now; its after side is drawn already.
  assert.equal(await history.keepRedrawnResult(), false, 'undone');
  await history.redo();
  assert.deepEqual(restores.at(-1), { drawn: false });
  // A checkpoint step's side is kept only while nothing is open.
  let release;
  const open = history.runUndoableCheckpoint('Open step', () => new Promise((resolve) => { release = resolve; }));
  await Promise.resolve();
  assert.equal(await history.keepRedrawnResult(), false, 'a step is open');
  release();
  await open;
  // A step without a caller's answer restores as before.
  await history.runUndoableCheckpoint('Reset', () => { rules = 4; });
  await history.undo();
  assert.deepEqual(restores.at(-1), { drawn: false });
  // N1 (review C): the before side is drawn as History captures it, not as
  // the caller asked for the step: a step that waits for an open one reads it
  // once that one settled.
  let owed = false;
  let settle;
  const waiting = history.runUndoable('Edit setting', () => new Promise((resolve) => { settle = resolve; }));
  await new Promise((resolve) => { setTimeout(resolve, 0); });
  assert.equal(typeof settle, 'function', 'the step is open');
  const reads = [];
  const recolored = history.runUndoableCheckpoint('Change feature color late', () => {
    rules = 5;
    figure = 'shown-5';
  }, { drawn: { before: () => { reads.push(owed); return !owed; }, change: true } });
  owed = true;
  rules = 6;
  settle();
  await waiting;
  await recolored;
  assert.deepEqual(reads, [true], 'read once, when the before side is captured');
  await history.undo();
  assert.deepEqual(restores.at(-1), { drawn: false }, 'owed when captured');
  await history.redo();
  assert.deepEqual(restores.at(-1), { drawn: false }, 'the after side follows the before side');
  console.log('history redrawn result tests passed');
}

// OV-872: what a shown rerender writes (the Legend rows Python drew) while an
// action's step is open, waiting for that rerender, is no change of a step:
// the action stays one step, its own writes stay in it, and its Undo keeps
// the drawn rows.
{
  let rules = 0;
  let legend = 'generated';
  let setting = 0;
  const history = createHistoryManager({
    buildIntent: () => ({ rules, legend, setting }),
    applyIntent: (intent) => { ({ rules, legend, setting } = intent); },
    buildCheckpoint: () => ({ rules, legend, setting }),
    applyCheckpoint: (checkpoint) => { ({ rules, legend, setting } = checkpoint); }
  });
  await history.captureBaseline('initial');
  /** @param {(shown: Promise<void>) => Promise<unknown>} fn */
  const whileDrawn = async (fn) => {
    let show = () => {};
    const shown = new Promise((resolve) => { show = resolve; });
    // The rerender was requested before the action's step began.
    const requestedAt = history.intentStepSerial();
    const action = history.runUndoable('Change feature color', () => fn(shown));
    await new Promise((resolve) => { setTimeout(resolve, 0); });
    history.adoptDrawnWrite(() => { legend = `drawn-${rules}`; }, { requestedAt });
    show();
    await action;
  };
  // An intent step: the action's write before and after the rerender.
  await whileDrawn(async (shown) => { setting = 1; await shown; rules = 1; });
  assert.equal(history.getUndoCount(), 1);
  await history.undo();
  assert.deepEqual({ rules, legend, setting }, { rules: 0, legend: 'drawn-0', setting: 0 });
  await history.redo();
  assert.deepEqual({ rules, legend, setting }, { rules: 1, legend: 'drawn-0', setting: 1 });
  // A rule commit: a checkpoint the action opens after the rerender.
  await whileDrawn(async (shown) => {
    await shown;
    await history.runUndoableCheckpoint('Change feature color', () => { rules = 2; });
  });
  assert.equal(history.getUndoCount(), 2, 'one step per action');
  await history.undo();
  assert.deepEqual({ rules, legend, setting }, { rules: 1, legend: 'drawn-1', setting: 1 });
  // With no step open, the write is in no step, as before.
  history.adoptDrawnWrite(() => { legend = 'drawn-again'; }, { requestedAt: history.intentStepSerial() });
  assert.equal(history.getUndoCount(), 1);
  console.log('history drawn write tests passed');
}

// D-40 (review S1 of #1091): the step that requested the rerender is still
// open when it binds (the label On path waits for the "Keep without label"
// choice). The drawn rows are its own: its Undo restores the rows from before
// its edit, as D-24 (OV-731) keeps a redraw's figure in the step that asked
// for it.
{
  let visible = false;
  let legend = ['generated'];
  const history = createHistoryManager({
    buildIntent: () => ({ visible, legend }),
    applyIntent: (intent) => { ({ visible, legend } = intent); },
    buildCheckpoint: () => ({ visible, legend }),
    applyCheckpoint: (checkpoint) => { ({ visible, legend } = checkpoint); }
  });
  await history.captureBaseline('initial');
  let show = () => {};
  const shown = new Promise((resolve) => { show = resolve; });
  /** @type {number | null} */
  let requestedAt = null;
  const action = history.runUndoable('Change label text', async () => {
    visible = true;
    requestedAt = history.intentStepSerial();
    await shown;
  });
  await new Promise((resolve) => { setTimeout(resolve, 0); });
  assert.notEqual(requestedAt, null, 'the step requested the rerender');
  history.adoptDrawnWrite(() => { legend = ['generated', 'shown feature']; }, { requestedAt: Number(requestedAt) });
  show();
  await action;
  assert.equal(history.getUndoCount(), 1);
  await history.undo();
  assert.deepEqual({ visible, legend }, { visible: false, legend: ['generated'] }, 'Undo restores the rows before the edit');
  await history.redo();
  assert.deepEqual({ visible, legend }, { visible: true, legend: ['generated', 'shown feature'] });

  // The limit (N1): arrays are adopted whole. A step that began after the
  // request and wrote the same array before the rerender binds loses that
  // write to the drawn array: its Undo does not revert it.
  const lateRequest = history.intentStepSerial();
  let release = () => {};
  const bound = new Promise((resolve) => { release = resolve; });
  const renaming = history.runUndoable('Rename Legend entry', async () => {
    legend = ['renamed', 'shown feature'];
    await bound;
  });
  await new Promise((resolve) => { setTimeout(resolve, 0); });
  history.adoptDrawnWrite(() => { legend = ['renamed', 'shown feature', 'drawn row']; }, { requestedAt: lateRequest });
  release();
  await renaming;
  assert.equal(history.getUndoCount(), 1, 'the rename is in no step: the drawn array replaced it');

  // N-A (review of #1091): a step that must first settle another owner's open
  // step is numbered when it begins, not after that wait, so a rerender its
  // own edit requests in the same task counts it as the requester.
  const pendingOther = await history.begin('Move legend', { owner: 'drag' });
  const opening = history.begin('Change setting', { owner: 'checkbox' });
  // The request is read in the task that began the step (a watcher flush).
  const ownRequest = history.intentStepSerial();
  const opened = await opening;
  assert.notEqual(opened, pendingOther, 'the other owner\'s step was settled');
  visible = false;
  history.adoptDrawnWrite(() => { legend = ['renamed', 'shown feature', 'own row']; }, { requestedAt: ownRequest });
  await history.commit(opened);
  await history.undo();
  assert.deepEqual({ visible, legend }, { visible: true, legend: ['renamed', 'shown feature', 'drawn row'] },
    'Undo restores the rows from before the step that requested them');
  console.log('history drawn write requester tests passed');
}

// D-07 (b) (Owner 2026-10-16 「確認を出す」): a change outside History (the
// Rich Feature Popup preference, D-03) makes the work unsaved without a step:
// Load Session asks before it replaces it. Save and Load clear the mark.
{
  let intent = { label: 'start' };
  const history = createHistoryManager({
    buildIntent: () => structuredClone(intent),
    applyIntent: (restored) => { intent = structuredClone(restored); },
    buildCheckpoint: () => structuredClone(intent),
    applyCheckpoint: (checkpoint) => { intent = structuredClone(checkpoint); }
  });
  await history.initializeIntentBaseline();
  await history.runUndoable('Change label text', () => { intent.label = 'edited'; });
  await history.undo();
  history.markSavePoint();
  assert.equal(await history.hasChangesSinceSavePoint(), false);
  history.markChangeOutsideHistory();
  assert.equal(await history.hasChangesSinceSavePoint(), true, 'a change outside History is unsaved');
  assert.deepEqual([history.getUndoCount(), history.getRedoCount()], [0, 1], 'it records no step and keeps Redo');
  history.markSavePoint();
  assert.equal(await history.hasChangesSinceSavePoint(), false, 'a Save clears the mark');
  history.markChangeOutsideHistory();
  await history.initializeIntentBaseline('Loaded session');
  assert.equal(await history.hasChangesSinceSavePoint(), false, 'a Load clears the mark');
  console.log('history change outside History tests passed');
}

// OV-1102: a bulk label edit replaces the saved label table (the rules a loaded
// Session sent) in its own step, so Undo returns both: the next Generate draws
// the saved labels and Save writes the table again. Redo replaces it again.
{
  const savedRow = {
    recordId: '*', featureType: '*', qualifier: 'product', valueRegex: '^tRNA-Phe$', labelText: 'SavedPhe'
  };
  const state = {
    form: {},
    adv: {},
    files: {},
    linearSeqs: [],
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] },
    results: ref([]),
    selectedResultIndex: ref(0),
    mode: ref('circular'),
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    downloadDpi: ref(300),
    canvasPadding: { top: 0, right: 0, bottom: 0, left: 0 },
    layoutPreferences: createLayoutPreferences(),
    extractedFeatures: ref([]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: {},
    featureVisibilityManualRules: [],
    featureOverrides: {},
    featureStrokeOverrides: {},
    canonicalLabelOverrideRows: ref([savedRow]),
    labelTextBulkOverrides: {},
    orthogroupDormantOverrides: { OG_1: { name: 'Saved group' } },
    orthogroups: ref([]),
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
    state: withDrawings(state),
    fileStore: createHistoryFileStore(),
    nextTick: async () => {},
    buildConfigData: () => ({ form: state.form, adv: state.adv })
  });
  const history = createHistoryManager({
    ...snapshots.drawingHistory('circular')
  });
  await history.initializeIntentBaseline();
  await history.runUndoable('Load label edits', () => {
    state.labelTextBulkOverrides.SavedPhe = 'MyLabel';
    state.canonicalLabelOverrideRows.value = [];
  });
  await history.undo();
  assert.deepEqual(
    [state.labelTextBulkOverrides, state.canonicalLabelOverrideRows.value],
    [{}, [savedRow]],
    'Undo returns the saved label table with the bulk label edits it replaced'
  );
  await history.redo();
  assert.deepEqual(
    [state.labelTextBulkOverrides, state.canonicalLabelOverrideRows.value],
    [{ SavedPhe: 'MyLabel' }, []]
  );
  // R15-21: Reset Settings removes the saved label table and the dormant
  // Similarity group names in its one step; Undo returns both.
  await history.undo();
  await history.runUndoableCheckpoint('Reset settings', () => {
    delete state.labelTextBulkOverrides.SavedPhe;
    state.canonicalLabelOverrideRows.value = [];
    delete state.orthogroupDormantOverrides.OG_1;
  });
  assert.deepEqual([state.canonicalLabelOverrideRows.value, state.orthogroupDormantOverrides], [[], {}]);
  await history.undo();
  assert.deepEqual(
    [state.canonicalLabelOverrideRows.value, state.orthogroupDormantOverrides],
    [[savedRow], { OG_1: { name: 'Saved group' } }],
    'Undo of Reset returns the saved label table and the dormant group names'
  );
  console.log('history saved label table tests passed');
}

// A3: one History manager per drawing behind a router. The managers share one
// byte budget, in which the steps of drawings that are not shown go first, and
// one file retention, which keeps a file while any manager names it (OV-112).
{
  const { createHistoryRouter } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));
  assert.equal(typeof createHistoryRouter, 'function');
  const drawingIds = ['circular', 'linear'];
  const fileA = makeFile('a.gb', 10);
  const fileB = makeFile('b.gb', 10);
  const makeRouter = (maxBytes) => {
    const fileStore = createHistoryFileStore();
    const shown = ref('linear');
    const values = { circular: '', linear: '' };
    const files = { circular: null, linear: null };
    const notices = [];
    const router = createHistoryRouter({
      drawingIds,
      activeDrawingId: () => shown.value,
      maxBytes,
      onHistoryLimit: (message) => { notices.push(message); },
      sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
      estimateFileBytes: fileStore.estimateBytes,
      retainOnlyFiles: fileStore.retainOnly,
      managerOptions: (drawingId) => {
        const buildState = () => ({
          value: values[drawingId],
          file: files[drawingId] ? fileStore.describeFile(files[drawingId]) : null
        });
        const applyState = (snapshot) => {
          values[drawingId] = snapshot.value;
          files[drawingId] = snapshot.file ? fileStore.restoreFile(snapshot.file) : null;
        };
        return { buildIntent: buildState, applyIntent: applyState, buildCheckpoint: buildState, applyCheckpoint: applyState };
      }
    });
    return { router, shown, values, files, fileStore, notices };
  };

  // The size of one step, measured with no limit.
  const probe = makeRouter(Number.MAX_SAFE_INTEGER);
  await probe.router.captureBaseline();
  await probe.router.runUndoable('Probe', () => { probe.values.linear = 'x'.repeat(400); });
  const stepBytes = probe.router.getDiagnostics().retainedEntryBytes;
  assert.ok(stepBytes > 0);

  // Two steps exceed the budget; one fits. The hidden drawing's step goes.
  const { router, shown, values, files, fileStore, notices } = makeRouter(Math.floor(stepBytes * 1.5));
  files.linear = fileA;
  await router.captureBaseline();
  // The file's id as History names it (`describeFile` gives a released file a new one).
  const fileAId = router.managerOf('linear').getCurrentIntent().file.fileId;
  await router.runUndoable('Linear edit', () => { values.linear = 'x'.repeat(400); files.linear = null; });
  assert.equal(router.canUndo(), true);
  shown.value = 'circular';
  assert.equal(router.canUndo(), false, 'the router answers for the shown drawing');
  await router.runUndoable('Circular edit', () => { values.circular = 'y'.repeat(400); });
  assert.equal(router.canUndo(), true, 'the shown drawing keeps its step');
  assert.equal(router.managerOf('linear').canUndo(), false, 'the hidden drawing loses its step first');
  assert.deepEqual(notices, ['Older undo history was discarded to stay within the history limit.']);
  assert.equal(fileStore.has(fileAId), false, 'no manager names file A any more');

  // A file stays while any manager names it: the Linear step names file B,
  // and a Circular step's release keeps it.
  shown.value = 'linear';
  files.linear = fileB;
  await router.captureBaseline();
  const fileBId = router.managerOf('linear').getCurrentIntent().file.fileId;
  await router.runUndoable('Remove B', () => { files.linear = null; });
  shown.value = 'circular';
  await router.runUndoable('Circular edit 2', () => { values.circular = 'z'; });
  assert.equal(fileStore.has(fileBId), true, 'the Linear step still names file B');
  shown.value = 'linear';
  await router.undo();
  assert.equal(files.linear?.name, 'b.gb');
  assert.equal(values.circular, 'z', 'Undo in Linear leaves Circular unchanged');

  // OV-1365: a pass that drops only the hidden drawing's Redo steps names the
  // redo history.
  const redoOnly = makeRouter(Math.floor(stepBytes * 1.5));
  await redoOnly.router.captureBaseline();
  await redoOnly.router.runUndoable('Linear edit', () => { redoOnly.values.linear = 'x'.repeat(400); });
  await redoOnly.router.undo();
  assert.equal(redoOnly.router.canRedo(), true);
  redoOnly.shown.value = 'circular';
  await redoOnly.router.runUndoable('Circular edit', () => { redoOnly.values.circular = 'y'.repeat(400); });
  assert.equal(redoOnly.router.managerOf('linear').canRedo(), false, 'the hidden drawing loses its Redo step');
  assert.deepEqual(redoOnly.notices, ['Redo history was discarded to stay within the history limit.']);

  // A3 review Low 4: the order follows the shown drawing, not the drawing
  // that commits. A step a hidden drawing records drops its own steps first.
  const hiddenCommit = makeRouter(Math.floor(stepBytes * 1.5));
  hiddenCommit.shown.value = 'circular';
  await hiddenCommit.router.captureBaseline();
  await hiddenCommit.router.runUndoable('Circular edit', () => { hiddenCommit.values.circular = 'y'.repeat(400); });
  await hiddenCommit.router.managerOf('linear').runUndoable('Hidden Linear edit', () => { hiddenCommit.values.linear = 'x'.repeat(400); });
  assert.equal(hiddenCommit.router.canUndo(), true, 'the shown drawing keeps its step');
  assert.equal(hiddenCommit.router.managerOf('linear').canUndo(), false, 'the hidden drawing that committed loses its step');
}

// OV-110 (D-15-49): the Result identities that History holds under a map keyed
// by identity (the composition offsets): the current intent's, an open step's,
// and those an Undo or Redo step names. A step the budget drops no longer
// holds its identities.
{
  const { createHistoryRouter } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));
  const shown = ref('circular');
  const ui = { circular: { offsets: {} }, linear: { offsets: {} } };
  const router = createHistoryRouter({
    drawingIds: ['circular', 'linear'],
    activeDrawingId: () => shown.value,
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    managerOptions: (drawingId) => {
      const build = () => ({ ui: JSON.parse(JSON.stringify(ui[drawingId])) });
      const apply = (snapshot) => { ui[drawingId] = JSON.parse(JSON.stringify(snapshot.ui)); };
      return { buildIntent: build, applyIntent: apply, buildCheckpoint: build, applyCheckpoint: apply, identityMapPath: ['ui', 'offsets'], maxActions: 2 };
    }
  });
  await router.captureBaseline();
  const held = () => [...router.heldIdentities()].sort();
  assert.deepEqual(held(), []);
  // R1 moved, then R1 replaced by R2 (no offsets yet): an Undo step still names R1.
  await router.runUndoable('Move R1', () => { ui.circular.offsets = { R1: [1, 2] }; });
  await router.runUndoable('Generate R2', () => { ui.circular.offsets = {}; });
  assert.deepEqual(held(), ['R1'], 'the Undo steps name R1');
  shown.value = 'linear';
  await router.runUndoable('Move L1', () => { ui.linear.offsets = { L1: [3, 4] }; });
  assert.deepEqual(held(), ['L1', 'R1'], 'every drawing\'s History');
  shown.value = 'circular';
  // Two more Circular steps drop the steps that named R1 (maxActions 2).
  await router.runUndoable('Move R2', () => { ui.circular.offsets = { R2: [5, 6] }; });
  await router.runUndoable('Move R2 again', () => { ui.circular.offsets = { R2: [7, 8] }; });
  assert.deepEqual(held(), ['L1', 'R2'], 'the dropped steps release R1');
  console.log('history held identity tests passed');
}

// OV-1339 (per-drawing History): the router's Undo or Redo is pending from its
// call, before the shown drawing settles its open step, to its end; the
// drawing switch reads `traversalPending` and waits for it.
{
  const { createHistoryRouter } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));
  const shown = ref('circular');
  const values = { circular: 0, linear: 0 };
  const router = createHistoryRouter({
    drawingIds: ['circular', 'linear'],
    activeDrawingId: () => shown.value,
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    managerOptions: (drawingId) => {
      const build = () => ({ value: values[drawingId] });
      const apply = (snapshot) => { values[drawingId] = snapshot.value; };
      return { buildIntent: build, applyIntent: apply, buildCheckpoint: build, applyCheckpoint: apply };
    }
  });
  await router.captureBaseline();
  await router.runUndoable('Increment', () => { values.circular = 1; });
  assert.equal(router.traversalPending(), false);
  const undoing = router.undo();
  assert.equal(router.traversalPending(), true, 'pending from the Undo call');
  await undoing;
  assert.equal(values.circular, 0);
  assert.equal(router.traversalPending(), false);
  const redoing = router.redo();
  assert.equal(router.traversalPending(), true, 'pending from the Redo call');
  await redoing;
  assert.equal(router.traversalPending(), false);
  console.log('history router traversal pending tests passed');
}
