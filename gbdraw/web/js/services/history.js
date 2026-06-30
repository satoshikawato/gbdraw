import { collectHistoryFileIds } from './history-files.js';

const DEFAULT_MAX_ACTIONS = 30;
const DEFAULT_MAX_BYTES = 200 * 1024 * 1024;

const makeBox = (value) => ({ value });

const estimateJsonBytes = (value) => {
  try {
    return JSON.stringify(value).length * 2;
  } catch (_err) {
    return DEFAULT_MAX_BYTES + 1;
  }
};

export const createHistoryManager = ({
  buildSnapshot,
  applySnapshot,
  snapshotSignature = (snapshot) => JSON.stringify(snapshot),
  fileStore = null,
  maxActions = DEFAULT_MAX_ACTIONS,
  maxBytes = DEFAULT_MAX_BYTES,
  makeRef = makeBox
} = {}) => {
  if (typeof buildSnapshot !== 'function') {
    throw new Error('createHistoryManager requires buildSnapshot.');
  }
  if (typeof applySnapshot !== 'function') {
    throw new Error('createHistoryManager requires applySnapshot.');
  }

  const undoStack = [];
  const redoStack = [];
  const revision = makeRef(0);
  const restoring = makeRef(false);
  const capturing = makeRef(false);
  const historyLimitMessage = makeRef('');
  let currentSnapshot = null;
  let currentSignature = '';
  let activeTransaction = null;

  const touch = () => {
    revision.value += 1;
  };

  const captureSnapshot = async () => {
    capturing.value = true;
    try {
      const snapshot = await buildSnapshot();
      const signature = snapshotSignature(snapshot);
      return { snapshot, signature };
    } finally {
      capturing.value = false;
    }
  };

  const collectReferencedFileIds = () => {
    const ids = new Set();
    if (currentSnapshot) collectHistoryFileIds(currentSnapshot, ids);
    undoStack.forEach((entry) => {
      collectHistoryFileIds(entry.before, ids);
      collectHistoryFileIds(entry.after, ids);
    });
    redoStack.forEach((entry) => {
      collectHistoryFileIds(entry.before, ids);
      collectHistoryFileIds(entry.after, ids);
    });
    return ids;
  };

  const estimateHistoryBytes = () => {
    const uniqueSnapshots = new Map();
    const addSnapshot = (snapshot) => {
      if (!snapshot) return;
      const signature = snapshotSignature(snapshot);
      if (!uniqueSnapshots.has(signature)) uniqueSnapshots.set(signature, snapshot);
    };

    addSnapshot(currentSnapshot);
    undoStack.forEach((entry) => {
      addSnapshot(entry.before);
      addSnapshot(entry.after);
    });
    redoStack.forEach((entry) => {
      addSnapshot(entry.before);
      addSnapshot(entry.after);
    });

    let bytes = 0;
    uniqueSnapshots.forEach((snapshot) => {
      bytes += estimateJsonBytes(snapshot);
    });
    if (fileStore?.estimateBytes) {
      bytes += fileStore.estimateBytes(collectReferencedFileIds());
    }
    return bytes;
  };

  const releaseUnreferencedFiles = () => {
    if (!fileStore?.retainOnly) return;
    fileStore.retainOnly(collectReferencedFileIds());
  };

  const enforceLimits = () => {
    let evicted = false;
    while (undoStack.length > maxActions) {
      undoStack.shift();
      evicted = true;
    }

    let bytes = estimateHistoryBytes();
    while (undoStack.length > 0 && bytes > maxBytes) {
      undoStack.shift();
      evicted = true;
      bytes = estimateHistoryBytes();
    }
    if (bytes > maxBytes) {
      redoStack.splice(0, redoStack.length);
      bytes = estimateHistoryBytes();
    }

    historyLimitMessage.value = evicted
      ? 'Older undo history was discarded to stay within the history limit.'
      : '';
    releaseUnreferencedFiles();
  };

  const setCurrent = (snapshot, signature = snapshotSignature(snapshot)) => {
    currentSnapshot = snapshot;
    currentSignature = signature;
  };

  const captureBaseline = async (_label = 'Baseline') => {
    if (restoring.value) return;
    const { snapshot, signature } = await captureSnapshot();
    undoStack.splice(0, undoStack.length);
    redoStack.splice(0, redoStack.length);
    setCurrent(snapshot, signature);
    enforceLimits();
    touch();
  };

  const begin = async (label = 'Edit', options = {}) => {
    if (restoring.value || capturing.value) return null;
    if (activeTransaction && !activeTransaction.closed) return activeTransaction;

    const { snapshot, signature } = await captureSnapshot();
    const tx = {
      label,
      before: snapshot,
      beforeSignature: signature,
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

  const commit = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    const { snapshot, signature } = await captureSnapshot();
    transaction.closed = true;
    if (activeTransaction === transaction) activeTransaction = null;

    if (signature === transaction.beforeSignature) {
      setCurrent(snapshot, signature);
      releaseUnreferencedFiles();
      touch();
      return false;
    }

    undoStack.push({
      label: transaction.label || options.label || 'Edit',
      before: transaction.before,
      beforeSignature: transaction.beforeSignature,
      after: snapshot,
      afterSignature: signature
    });
    redoStack.splice(0, redoStack.length);
    setCurrent(snapshot, signature);
    enforceLimits();
    touch();
    return true;
  };

  const runUndoable = async (label, fn, options = {}) => {
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn();

    const usesActiveTransaction = Boolean(activeTransaction && !activeTransaction.closed);
    const tx = usesActiveTransaction
      ? activeTransaction
      : await begin(label, options);
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

  const applySnapshotWithFlag = async (snapshot) => {
    restoring.value = true;
    try {
      await applySnapshot(snapshot);
    } finally {
      restoring.value = false;
    }
  };

  const undo = async () => {
    if (restoring.value || undoStack.length === 0) return false;
    const entry = undoStack.pop();
    redoStack.push(entry);
    await applySnapshotWithFlag(entry.before);
    setCurrent(entry.before, entry.beforeSignature);
    enforceLimits();
    touch();
    return true;
  };

  const redo = async () => {
    if (restoring.value || redoStack.length === 0) return false;
    const entry = redoStack.pop();
    undoStack.push(entry);
    await applySnapshotWithFlag(entry.after);
    setCurrent(entry.after, entry.afterSignature);
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
    cancel,
    captureBaseline,
    canRedo,
    canUndo,
    commit,
    getCurrentSnapshot: () => currentSnapshot,
    getCurrentSignature: () => currentSignature,
    getRedoCount: () => redoStack.length,
    getUndoCount: () => undoStack.length,
    historyLimitMessage,
    redo,
    redoLabel,
    restoring,
    capturing,
    revision,
    runUndoable,
    undo,
    undoLabel
  };
};
