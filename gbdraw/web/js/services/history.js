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
  signatureFor = (value) => JSON.stringify(value),
  fileStore = null,
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
    byteEstimateComputations: 0,
    historySvgBytes: 0
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
        fileIds: collectHistoryFileIds(intent, new Set())
      };
    } finally {
      capturing.value = false;
    }
  };

  const finishCheckpointCapture = (checkpoint) => {
    const signature = computeSignature(checkpoint, 'artifact-checkpoint');
    const svgBytes = checkpointSvgBytes(checkpoint);
    diagnostics.historySvgBytes += svgBytes;
    emitHistoryDiagnostic({ type: 'capture', scope: 'artifact-checkpoint', svgBytes });
    return {
      checkpoint,
      signature,
      byteSize: estimateSignedBytes(signature, 'artifact-checkpoint'),
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
      emitHistoryDiagnostic({ type: 'begin', scope: 'artifact-checkpoint', label });
      return {
        label,
        before,
        closed: false,
        source: options.source || ''
      };
    };
    const captureBefore = () => {
      const before = captureCheckpoint();
      return isPromiseLike(before)
        ? Promise.resolve(before).then(createTransaction)
        : createTransaction(before);
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
    return true;
  };

  const runUndoableCheckpoint = (label, fn, options = {}) => {
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn();
    const execute = async (tx) => {
      try {
        const result = await fn();
        await commitCheckpoint(tx, options);
        return result;
      } catch (error) {
        if (tx) tx.closed = true;
        throw error;
      }
    };
    const transaction = beginCheckpoint(label, options);
    return isPromiseLike(transaction)
      ? Promise.resolve(transaction).then(execute)
      : execute(transaction);
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
      fileIds: collectHistoryFileIds(nextIntent, new Set())
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
    beginCheckpoint,
    cancel,
    captureBaseline,
    canRedo,
    canUndo,
    commit,
    commitCheckpoint,
    getCurrentCheckpoint: () => currentCheckpoint,
    getCurrentCheckpointSignature: () => currentCheckpointSignature,
    getCurrentIntent: () => currentIntent,
    getCurrentIntentSignature: () => currentIntentSignature,
    getDiagnostics: () => ({ ...diagnostics, retainedEntryBytes: totalEntryBytes }),
    getRedoCount: () => redoStack.length,
    getUndoCount: () => undoStack.length,
    historyLimitMessage,
    redo,
    redoLabel,
    restoring,
    capturing,
    revision,
    runUndoable,
    runUndoableCheckpoint,
    runUndoableCommand,
    undo,
    undoLabel
  };
};
