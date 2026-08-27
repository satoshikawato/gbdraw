import { collectHistoryFileIds } from './history-files.js';

const DEFAULT_MAX_ACTIONS = 30;
const DEFAULT_MAX_BYTES = 200 * 1024 * 1024;

const makeBox = (value) => ({ value });
const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);
const isPromiseLike = (value) => (
  Boolean(value) && typeof value.then === 'function'
);

const cloneIntentValue = (value) => {
  if (Array.isArray(value)) return value.map((entry) => cloneIntentValue(entry));
  if (!isPlainObject(value)) return value;
  return Object.fromEntries(
    Object.entries(value).map(([key, entry]) => [key, cloneIntentValue(entry)])
  );
};

const sameJsonValue = (left, right) => {
  if (Object.is(left, right)) return true;
  if (typeof left !== typeof right) return false;
  if (left === null || right === null) return false;
  if (typeof left !== 'object') return false;
  return JSON.stringify(left) === JSON.stringify(right);
};

const buildIntentPatch = (before, after, path = [], changes = []) => {
  if (Object.is(before, after)) return changes;

  const beforeObject = isPlainObject(before);
  const afterObject = isPlainObject(after);
  if (beforeObject && afterObject) {
    const keys = new Set([...Object.keys(before), ...Object.keys(after)]);
    keys.forEach((key) => {
      const beforeHas = hasOwn(before, key);
      const afterHas = hasOwn(after, key);
      if (!beforeHas || !afterHas) {
        changes.push({
          path: [...path, key],
          before: beforeHas ? cloneIntentValue(before[key]) : undefined,
          after: afterHas ? cloneIntentValue(after[key]) : undefined,
          beforeHas,
          afterHas
        });
        return;
      }
      buildIntentPatch(before[key], after[key], [...path, key], changes);
    });
    return changes;
  }

  if (!sameJsonValue(before, after)) {
    changes.push({
      path,
      before: cloneIntentValue(before),
      after: cloneIntentValue(after),
      beforeHas: true,
      afterHas: true
    });
  }
  return changes;
};

const applyIntentPatch = (source, changes, direction) => {
  const target = cloneIntentValue(source || {});
  const valueKey = direction === 'undo' ? 'before' : 'after';
  const presentKey = direction === 'undo' ? 'beforeHas' : 'afterHas';

  changes.forEach((change) => {
    const path = Array.isArray(change.path) ? change.path : [];
    if (path.length === 0) return;
    let owner = target;
    for (let index = 0; index < path.length - 1; index += 1) {
      const key = path[index];
      if (!isPlainObject(owner[key]) && !Array.isArray(owner[key])) owner[key] = {};
      owner = owner[key];
    }
    const key = path[path.length - 1];
    if (!change[presentKey]) {
      delete owner[key];
    } else {
      owner[key] = cloneIntentValue(change[valueKey]);
    }
  });

  return target;
};

const emitHistoryDiagnostic = (event) => {
  const callback = globalThis.__GBDRAW_TEST_HOOKS__?.onHistoryDiagnostic;
  if (typeof callback !== 'function') return;
  try {
    callback(event);
  } catch (_error) {
    // Test instrumentation is observational and must not own application behavior.
  }
};

const checkpointSvgBytes = (checkpoint) => (
  Array.isArray(checkpoint?.results)
    ? checkpoint.results.reduce(
        (total, result) => total + String(result?.content || '').length * 2,
        0
      )
    : 0
);

export const createHistoryManager = ({
  buildIntent,
  applyIntent,
  buildCheckpoint,
  applyCheckpoint,
  captureGeneratedArtifactHandle = null,
  restoreGeneratedArtifactHandle = null,
  compareGeneratedArtifactHandles = null,
  signatureFor = (value) => JSON.stringify(value),
  fileStore = null,
  collectCurrentFileIds = null,
  maxActions = DEFAULT_MAX_ACTIONS,
  maxBytes = DEFAULT_MAX_BYTES,
  makeRef = makeBox
} = {}) => {
  if (typeof buildIntent !== 'function') {
    throw new Error('createHistoryManager requires buildIntent.');
  }
  if (typeof applyIntent !== 'function') {
    throw new Error('createHistoryManager requires applyIntent.');
  }
  if (typeof buildCheckpoint !== 'function') {
    throw new Error('createHistoryManager requires buildCheckpoint.');
  }
  if (typeof applyCheckpoint !== 'function') {
    throw new Error('createHistoryManager requires applyCheckpoint.');
  }

  const undoStack = [];
  const redoStack = [];
  const revision = makeRef(0);
  const restoring = makeRef(false);
  const capturing = makeRef(false);
  const historyLimitMessage = makeRef('');
  const diagnostics = {
    intentBuilds: 0,
    artifactCheckpointBuilds: 0,
    signatureComputations: 0,
    intentSignatureComputations: 0,
    artifactCheckpointSignatureComputations: 0,
    byteEstimateComputations: 0,
    historySvgBytes: 0,
    checkpointEstimatedBytes: 0,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0,
    manualCancelFullArtifactSnapshotBuildCount: 0,
    artifactHandleBeforeBuildCount: 0,
    artifactHandleAfterBuildCount: 0,
    artifactFingerprintComparisonCount: 0,
    artifactReplacementHistoryEntryCount: 0
  };
  let currentIntent = null;
  let currentIntentSignature = '';
  let currentCheckpoint = null;
  let currentCheckpointSignature = '';
  let currentFileIds = new Set();
  let totalEntryBytes = 0;
  let activeTransaction = null;

  const touch = () => {
    revision.value += 1;
  };

  const computeSignature = (value, scope) => {
    diagnostics.signatureComputations += 1;
    if (scope === 'artifact-checkpoint') {
      diagnostics.artifactCheckpointSignatureComputations += 1;
    } else if (String(scope || '').startsWith('intent')) {
      diagnostics.intentSignatureComputations += 1;
    }
    const signature = signatureFor(value);
    emitHistoryDiagnostic({ type: 'signature', scope, bytes: String(signature || '').length * 2 });
    return signature;
  };

  const estimateSignedBytes = (signature, scope) => {
    diagnostics.byteEstimateComputations += 1;
    const bytes = String(signature || '').length * 2;
    emitHistoryDiagnostic({ type: 'size', scope, bytes });
    return bytes;
  };

  const notifyCheckpointCapture = (options, phase) => {
    const callback = options?.onCheckpointCapture;
    if (typeof callback !== 'function') return;
    try {
      callback({ phase, diagnostics: { ...diagnostics } });
    } catch (_error) {
      // Diagnostic observers cannot own History behavior.
    }
  };

  const collectIntentFileIds = (intent) => {
    const fileIds = collectHistoryFileIds(intent, new Set());
    if (typeof collectCurrentFileIds !== 'function') return fileIds;
    const currentIds = collectCurrentFileIds();
    if (currentIds && typeof currentIds[Symbol.iterator] === 'function') {
      for (const fileId of currentIds) fileIds.add(fileId);
    }
    return fileIds;
  };

  const captureIntent = async () => {
    capturing.value = true;
    diagnostics.intentBuilds += 1;
    try {
      const intent = await buildIntent();
      const signature = computeSignature(intent, 'intent');
      emitHistoryDiagnostic({ type: 'capture', scope: 'intent' });
      return {
        intent,
        signature,
        fileIds: collectIntentFileIds(intent)
      };
    } finally {
      capturing.value = false;
    }
  };

  const finishCheckpointCapture = (checkpoint) => {
    const signature = computeSignature(checkpoint, 'artifact-checkpoint');
    const svgBytes = checkpointSvgBytes(checkpoint);
    const byteSize = estimateSignedBytes(signature, 'artifact-checkpoint');
    diagnostics.historySvgBytes += svgBytes;
    diagnostics.checkpointEstimatedBytes += byteSize;
    emitHistoryDiagnostic({ type: 'capture', scope: 'artifact-checkpoint', svgBytes });
    return {
      checkpoint,
      signature,
      byteSize,
      fileIds: collectHistoryFileIds(checkpoint, new Set())
    };
  };

  const captureCheckpoint = () => {
    capturing.value = true;
    diagnostics.artifactCheckpointBuilds += 1;
    try {
      const checkpoint = buildCheckpoint();
      if (isPromiseLike(checkpoint)) {
        return Promise.resolve(checkpoint)
          .then(finishCheckpointCapture)
          .finally(() => {
            capturing.value = false;
          });
      }
      const record = finishCheckpointCapture(checkpoint);
      capturing.value = false;
      return record;
    } catch (error) {
      capturing.value = false;
      throw error;
    }
  };

  const setCurrentIntent = ({ intent, signature, fileIds }) => {
    currentIntent = intent;
    currentIntentSignature = signature;
    currentFileIds = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
  };

  const setCurrentCheckpoint = ({ checkpoint, signature }) => {
    currentCheckpoint = checkpoint;
    currentCheckpointSignature = signature;
  };

  const clearCurrentCheckpoint = () => {
    currentCheckpoint = null;
    currentCheckpointSignature = '';
  };

  const entryFileIds = (entry) => (
    entry?.fileIds instanceof Set ? entry.fileIds : new Set()
  );

  const collectRetainedFileIds = () => {
    const ids = new Set(currentFileIds);
    undoStack.forEach((entry) => entryFileIds(entry).forEach((id) => ids.add(id)));
    redoStack.forEach((entry) => entryFileIds(entry).forEach((id) => ids.add(id)));
    return ids;
  };

  const retainedFileBytes = () => (
    fileStore?.estimateBytes ? fileStore.estimateBytes(collectRetainedFileIds()) : 0
  );

  const releaseUnreferencedFiles = () => {
    if (!fileStore?.retainOnly) return;
    fileStore.retainOnly(collectRetainedFileIds());
  };

  const addEntryBytes = (entry) => {
    totalEntryBytes += Number(entry?.byteSize) || 0;
  };

  const removeEntryBytes = (entry) => {
    totalEntryBytes = Math.max(0, totalEntryBytes - (Number(entry?.byteSize) || 0));
  };

  const clearStack = (stack) => {
    stack.forEach(removeEntryBytes);
    stack.splice(0, stack.length);
  };

  const pushUndoEntry = (entry) => {
    undoStack.push(entry);
    addEntryBytes(entry);
  };

  const clearRedo = () => clearStack(redoStack);

  const enforceLimits = () => {
    let evicted = false;
    while (undoStack.length > maxActions) {
      removeEntryBytes(undoStack.shift());
      evicted = true;
    }

    while (undoStack.length > 0 && totalEntryBytes + retainedFileBytes() > maxBytes) {
      removeEntryBytes(undoStack.shift());
      evicted = true;
    }
    if (totalEntryBytes + retainedFileBytes() > maxBytes) clearRedo();

    historyLimitMessage.value = evicted
      ? 'Older undo history was discarded to stay within the history limit.'
      : '';
    releaseUnreferencedFiles();
  };

  const captureBaseline = async (_label = 'Baseline') => {
    if (restoring.value) return;
    const checkpointRecord = await captureCheckpoint();
    const intentRecord = await captureIntent();
    clearStack(undoStack);
    clearRedo();
    setCurrentCheckpoint(checkpointRecord);
    setCurrentIntent(intentRecord);
    enforceLimits();
    touch();
  };

  const begin = async (label = 'Edit', options = {}) => {
    if (restoring.value || capturing.value) return null;
    if (activeTransaction && !activeTransaction.closed) return activeTransaction;

    const before = await captureIntent();
    emitHistoryDiagnostic({ type: 'begin', scope: 'intent', label });
    const tx = {
      label,
      before: before.intent,
      beforeSignature: before.signature,
      beforeFileIds: before.fileIds,
      closed: false,
      source: options.source || ''
    };
    activeTransaction = tx;
    return tx;
  };

  const cancel = (transaction) => {
    if (!transaction) return;
    transaction.closed = true;
    if (activeTransaction === transaction) activeTransaction = null;
  };

  const initializeIntentBaseline = async (_label = 'Intent baseline') => {
    if (restoring.value) return false;
    if (activeTransaction && !activeTransaction.closed) cancel(activeTransaction);
    const intentRecord = await captureIntent();
    clearStack(undoStack);
    clearRedo();
    currentCheckpoint = null;
    currentCheckpointSignature = '';
    setCurrentIntent(intentRecord);
    enforceLimits();
    touch();
    return true;
  };

  const buildPatchEntry = (transaction, afterRecord, options = {}) => {
    if (afterRecord.signature === transaction.beforeSignature) return null;
    const changes = buildIntentPatch(transaction.before, afterRecord.intent);
    if (changes.length === 0) return null;
    const metadataSignature = JSON.stringify(changes);
    diagnostics.byteEstimateComputations += 1;
    emitHistoryDiagnostic({ type: 'size', scope: 'intent-entry', bytes: metadataSignature.length * 2 });
    const fileIds = new Set(transaction.beforeFileIds || []);
    afterRecord.fileIds.forEach((id) => fileIds.add(id));
    return {
      type: 'intent',
      label: transaction.label || options.label || 'Edit',
      changes,
      byteSize: metadataSignature.length * 2,
      fileIds
    };
  };

  const commit = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    const afterRecord = await captureIntent();
    transaction.closed = true;
    if (activeTransaction === transaction) activeTransaction = null;
    setCurrentIntent(afterRecord);

    const entry = buildPatchEntry(transaction, afterRecord, options);
    if (!entry) {
      emitHistoryDiagnostic({ type: 'commit', scope: 'intent', label: transaction.label, created: false });
      releaseUnreferencedFiles();
      touch();
      return false;
    }

    pushUndoEntry(entry);
    emitHistoryDiagnostic({
      type: 'commit',
      scope: 'intent',
      label: entry.label,
      created: true,
      changeCount: entry.changes.length,
      bytes: entry.byteSize
    });
    clearRedo();
    enforceLimits();
    touch();
    return true;
  };

  const runUndoable = async (label, fn, options = {}) => {
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn();

    const usesActiveTransaction = Boolean(activeTransaction && !activeTransaction.closed);
    const tx = usesActiveTransaction ? activeTransaction : await begin(label, options);
    if (usesActiveTransaction && tx) tx.deferAdapterCommit = true;

    try {
      const result = await fn();
      await commit(tx, options);
      return result;
    } catch (error) {
      cancel(tx);
      throw error;
    }
  };

  const beginCheckpoint = (label = 'Artifact change', options = {}) => {
    if (restoring.value || capturing.value) return null;
    const createTransaction = (before) => {
      if (currentCheckpoint !== null) setCurrentCheckpoint(before);
      emitHistoryDiagnostic({ type: 'begin', scope: 'artifact-checkpoint', label });
      return {
        label,
        before,
        closed: false,
        source: options.source || ''
      };
    };
    const captureBefore = () => {
      notifyCheckpointCapture(options, 'before-start');
      const before = captureCheckpoint();
      return isPromiseLike(before)
        ? Promise.resolve(before).then((record) => {
            const transaction = createTransaction(record);
            notifyCheckpointCapture(options, 'before-end');
            return transaction;
          })
        : (() => {
            const transaction = createTransaction(before);
            notifyCheckpointCapture(options, 'before-end');
            return transaction;
          })();
    };

    if (!activeTransaction || activeTransaction.closed) return captureBefore();

    const pendingIntent = activeTransaction;
    const previousDeferred = Boolean(pendingIntent.deferAdapterCommit);
    pendingIntent.deferAdapterCommit = true;
    return Promise.resolve(commit(pendingIntent))
      .then(captureBefore)
      .catch((error) => {
        if (!pendingIntent.closed) pendingIntent.deferAdapterCommit = previousDeferred;
        throw error;
      });
  };

  const commitCheckpoint = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    notifyCheckpointCapture(options, 'after-start');
    const after = await captureCheckpoint();
    const intent = await captureIntent();
    transaction.closed = true;
    setCurrentCheckpoint(after);
    setCurrentIntent(intent);

    if (after.signature === transaction.before.signature) {
      emitHistoryDiagnostic({
        type: 'commit',
        scope: 'artifact-checkpoint',
        label: transaction.label,
        created: false
      });
      releaseUnreferencedFiles();
      touch();
      notifyCheckpointCapture(options, 'after-end');
      return false;
    }

    const fileIds = new Set(transaction.before.fileIds);
    after.fileIds.forEach((id) => fileIds.add(id));
    const entry = {
      type: 'checkpoint',
      label: transaction.label || options.label || 'Artifact change',
      before: transaction.before,
      after,
      byteSize: transaction.before.byteSize + after.byteSize,
      fileIds
    };
    pushUndoEntry(entry);
    emitHistoryDiagnostic({
      type: 'commit',
      scope: 'artifact-checkpoint',
      label: entry.label,
      created: true,
      bytes: entry.byteSize
    });
    clearRedo();
    enforceLimits();
    touch();
    notifyCheckpointCapture(options, 'after-end');
    return true;
  };

  const runUndoableCheckpoint = (label, fn, options = {}) => {
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn();
    const execute = async (tx) => {
      try {
        const result = await fn();
        if (
          typeof options.shouldCommit === 'function'
          && !options.shouldCommit(result)
        ) {
          if (tx) tx.closed = true;
          releaseUnreferencedFiles();
          return result;
        }
        await commitCheckpoint(tx, options);
        return result;
      } catch (error) {
        if (tx) tx.closed = true;
        releaseUnreferencedFiles();
        throw error;
      }
    };
    const transaction = beginCheckpoint(label, options);
    return isPromiseLike(transaction)
      ? Promise.resolve(transaction).then(execute)
      : execute(transaction);
  };

  const finishArtifactHandleCapture = (handle, phase) => {
    if (!handle || typeof handle !== 'object') {
      throw new Error('Generated artifact capture did not return a handle.');
    }
    const retainedBytes = Number(handle.retainedBytes);
    if (!Number.isFinite(retainedBytes) || retainedBytes < 0) {
      throw new Error('Generated artifact handle has an invalid retained-byte estimate.');
    }
    emitHistoryDiagnostic({
      type: 'capture',
      scope: 'artifact-replacement',
      phase,
      bytes: retainedBytes
    });
    return handle;
  };

  const captureArtifactHandle = (phase) => {
    if (typeof captureGeneratedArtifactHandle !== 'function') {
      throw new Error('Generated artifact replacement requires a capture owner.');
    }
    capturing.value = true;
    if (phase === 'before') diagnostics.artifactHandleBeforeBuildCount += 1;
    else diagnostics.artifactHandleAfterBuildCount += 1;
    try {
      const handle = captureGeneratedArtifactHandle({ phase });
      if (isPromiseLike(handle)) {
        return Promise.resolve(handle)
          .then((captured) => finishArtifactHandleCapture(captured, phase))
          .finally(() => {
            capturing.value = false;
          });
      }
      const captured = finishArtifactHandleCapture(handle, phase);
      capturing.value = false;
      return captured;
    } catch (error) {
      capturing.value = false;
      throw error;
    }
  };

  const beginArtifactReplacement = async (label = 'Generate diagram', options = {}) => {
    if (restoring.value || capturing.value) return null;
    if (activeTransaction && !activeTransaction.closed) {
      const pendingIntent = activeTransaction;
      const previousDeferred = Boolean(pendingIntent.deferAdapterCommit);
      pendingIntent.deferAdapterCommit = true;
      try {
        await commit(pendingIntent);
      } catch (error) {
        if (!pendingIntent.closed) pendingIntent.deferAdapterCommit = previousDeferred;
        throw error;
      }
    }
    notifyCheckpointCapture(options, 'before-start');
    const before = await captureArtifactHandle('before');
    notifyCheckpointCapture(options, 'before-end');
    emitHistoryDiagnostic({ type: 'begin', scope: 'artifact-replacement', label });
    return {
      label,
      before,
      closed: false,
      source: options.source || ''
    };
  };

  const sameArtifactHandles = (before, after) => {
    diagnostics.artifactFingerprintComparisonCount += 1;
    const same = typeof compareGeneratedArtifactHandles === 'function'
      ? Boolean(compareGeneratedArtifactHandles(before, after))
      : Boolean(
          before?.identity?.fingerprint
          && before.identity.fingerprint === after?.identity?.fingerprint
          && before.identity.compactSignature === after?.identity?.compactSignature
        );
    emitHistoryDiagnostic({
      type: 'compare',
      scope: 'artifact-replacement',
      equal: same
    });
    return same;
  };

  const artifactHandleFileIds = (handle) => new Set(handle?.fileIds || []);

  const commitAppliedArtifactReplacement = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    notifyCheckpointCapture(options, 'after-start');
    const after = await captureArtifactHandle('after');
    const intent = await captureIntent();
    transaction.closed = true;
    setCurrentIntent(intent);
    clearCurrentCheckpoint();

    if (sameArtifactHandles(transaction.before, after)) {
      emitHistoryDiagnostic({
        type: 'commit',
        scope: 'artifact-replacement',
        label: transaction.label,
        created: false
      });
      releaseUnreferencedFiles();
      touch();
      notifyCheckpointCapture(options, 'after-end');
      return false;
    }

    const fileIds = artifactHandleFileIds(transaction.before);
    artifactHandleFileIds(after).forEach((id) => fileIds.add(id));
    const entry = {
      type: 'artifact-replacement',
      label: transaction.label || options.label || 'Generate diagram',
      before: transaction.before,
      after,
      byteSize: Math.max(
        1,
        (Number(transaction.before.retainedBytes) || 0) + (Number(after.retainedBytes) || 0)
      ),
      fileIds
    };
    diagnostics.byteEstimateComputations += 1;
    diagnostics.artifactReplacementHistoryEntryCount += 1;
    emitHistoryDiagnostic({
      type: 'size',
      scope: 'artifact-replacement',
      bytes: entry.byteSize
    });
    pushUndoEntry(entry);
    clearRedo();
    enforceLimits();
    touch();
    emitHistoryDiagnostic({
      type: 'commit',
      scope: 'artifact-replacement',
      label: entry.label,
      created: true,
      bytes: entry.byteSize
    });
    notifyCheckpointCapture(options, 'after-end');
    return true;
  };

  const restoreArtifactHandle = async (handle) => {
    if (typeof restoreGeneratedArtifactHandle !== 'function') {
      throw new Error('Generated artifact replacement requires a restore owner.');
    }
    restoring.value = true;
    try {
      await restoreGeneratedArtifactHandle(handle, {
        clearFailedGeneratePresentation: true
      });
    } finally {
      restoring.value = false;
    }
    clearCurrentCheckpoint();
    await refreshCurrentIntent();
  };

  const runUndoableArtifactReplacement = async (label, fn, options = {}) => {
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn(null);
    const transaction = await beginArtifactReplacement(label, options);
    try {
      const result = await fn(transaction?.before || null);
      if (
        typeof options.shouldCommit === 'function'
        && !options.shouldCommit(result)
      ) {
        if (transaction) transaction.closed = true;
        await refreshCurrentIntent();
        releaseUnreferencedFiles();
        touch();
        return result;
      }
      await commitAppliedArtifactReplacement(transaction, options);
      return result;
    } catch (error) {
      if (transaction) {
        transaction.closed = true;
        await restoreArtifactHandle(transaction.before);
      }
      releaseUnreferencedFiles();
      throw error;
    }
  };

  const normalizeCommand = (label, command) => {
    if (!command || typeof command !== 'object' || command.noop) return null;
    if (typeof command.apply !== 'function') {
      throw new Error('Undoable command requires apply().');
    }
    if (typeof command.revert !== 'function') {
      throw new Error('Undoable command requires revert().');
    }
    diagnostics.byteEstimateComputations += 1;
    const byteSize = typeof command.estimateBytes === 'function'
      ? Number(command.estimateBytes()) || 0
      : JSON.stringify(command.metadata || {}).length * 2;
    emitHistoryDiagnostic({ type: 'size', scope: 'command', bytes: byteSize });
    return {
      type: 'command',
      label: command.label || label || 'Edit',
      apply: command.apply,
      revert: command.revert,
      byteSize,
      fileIds: collectHistoryFileIds(command.metadata || {}, new Set())
    };
  };

  const commandSucceeded = (result) => result !== false;

  const warnCommandNotApplied = (entry, action) => {
    const label = entry?.label ? ` "${entry.label}"` : '';
    console.warn(`Undo/redo command${label} could not ${action}; history was left unchanged.`);
  };

  const refreshCurrentIntent = async () => {
    setCurrentIntent(await captureIntent());
  };

  const runUndoableCommand = async (label, buildCommand) => {
    if (typeof buildCommand !== 'function') return false;
    if (!restoring.value && !capturing.value && activeTransaction && !activeTransaction.closed) {
      await commit(activeTransaction);
    }
    const command = normalizeCommand(label, await buildCommand());
    if (!command) return false;

    if (restoring.value || capturing.value) {
      const applied = commandSucceeded(await command.apply());
      if (applied) touch();
      return applied;
    }

    const applied = commandSucceeded(await command.apply());
    if (!applied) {
      warnCommandNotApplied(command, 'apply');
      return false;
    }
    await refreshCurrentIntent();
    pushUndoEntry(command);
    clearRedo();
    enforceLimits();
    touch();
    return true;
  };

  const applyIntentEntry = async (entry, direction) => {
    const nextIntent = applyIntentPatch(currentIntent || {}, entry.changes, direction);
    restoring.value = true;
    try {
      await applyIntent(nextIntent, { changes: entry.changes, direction });
    } finally {
      restoring.value = false;
    }
    const signature = computeSignature(nextIntent, 'intent-restore');
    setCurrentIntent({
      intent: nextIntent,
      signature,
      fileIds: collectIntentFileIds(nextIntent)
    });
  };

  const applyCheckpointEntry = async (record) => {
    restoring.value = true;
    try {
      await applyCheckpoint(record.checkpoint);
    } finally {
      restoring.value = false;
    }
    setCurrentCheckpoint(record);
    await refreshCurrentIntent();
  };

  const applyArtifactReplacementEntry = async (handle) => {
    await restoreArtifactHandle(handle);
  };

  const applyCommandWithFlag = async (entry, direction) => {
    restoring.value = true;
    try {
      const result = direction === 'undo' ? await entry.revert() : await entry.apply();
      return commandSucceeded(result);
    } finally {
      restoring.value = false;
    }
  };

  const undo = async () => {
    if (restoring.value || undoStack.length === 0) return false;
    const entry = undoStack[undoStack.length - 1];
    if (entry.type === 'command') {
      const reverted = await applyCommandWithFlag(entry, 'undo');
      if (!reverted) {
        warnCommandNotApplied(entry, 'undo');
        return false;
      }
      await refreshCurrentIntent();
    } else if (entry.type === 'checkpoint') {
      await applyCheckpointEntry(entry.before);
    } else if (entry.type === 'artifact-replacement') {
      await applyArtifactReplacementEntry(entry.before);
    } else {
      await applyIntentEntry(entry, 'undo');
    }
    undoStack.pop();
    redoStack.push(entry);
    enforceLimits();
    touch();
    return true;
  };

  const redo = async () => {
    if (restoring.value || redoStack.length === 0) return false;
    const entry = redoStack[redoStack.length - 1];
    if (entry.type === 'command') {
      const applied = await applyCommandWithFlag(entry, 'redo');
      if (!applied) {
        warnCommandNotApplied(entry, 'redo');
        return false;
      }
      await refreshCurrentIntent();
    } else if (entry.type === 'checkpoint') {
      await applyCheckpointEntry(entry.after);
    } else if (entry.type === 'artifact-replacement') {
      await applyArtifactReplacementEntry(entry.after);
    } else {
      await applyIntentEntry(entry, 'redo');
    }
    redoStack.pop();
    undoStack.push(entry);
    enforceLimits();
    touch();
    return true;
  };

  const canUndo = () => undoStack.length > 0 && !restoring.value;
  const canRedo = () => redoStack.length > 0 && !restoring.value;
  const undoLabel = () => (canUndo() ? undoStack[undoStack.length - 1].label : '');
  const redoLabel = () => (canRedo() ? redoStack[redoStack.length - 1].label : '');

  return {
    begin,
    beginArtifactReplacement,
    beginCheckpoint,
    cancel,
    captureBaseline,
    canRedo,
    canUndo,
    commit,
    commitAppliedArtifactReplacement,
    commitCheckpoint,
    getCurrentCheckpoint: () => currentCheckpoint,
    getCurrentCheckpointSignature: () => currentCheckpointSignature,
    getCurrentIntent: () => currentIntent,
    getCurrentIntentSignature: () => currentIntentSignature,
    getDiagnostics: () => ({ ...diagnostics, retainedEntryBytes: totalEntryBytes }),
    getRedoCount: () => redoStack.length,
    getUndoCount: () => undoStack.length,
    historyLimitMessage,
    initializeIntentBaseline,
    redo,
    redoLabel,
    restoring,
    capturing,
    revision,
    runUndoable,
    runUndoableArtifactReplacement,
    runUndoableCheckpoint,
    runUndoableCommand,
    undo,
    undoLabel
  };
};
