// @ts-check
import { collectHistoryFileIds } from './history-files.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

const DEFAULT_MAX_ACTIONS = 30;
const DEFAULT_MAX_BYTES = 200 * 1024 * 1024;
// The Session availability wording for an edit still being applied.
const HISTORY_EDIT_BUSY = Object.freeze({
  status: 'busy',
  reason: 'Applying an edit. Retry after the edit finishes.'
});

/**
 * @template T
 * @param {T} value
 * @returns {{ value: T }}
 */
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

/**
 * The artifact a generated-artifact handle names: its fingerprint and compact
 * signature, which the default handle comparison also reads. A handle without
 * a fingerprint names only itself.
 * @param {{ identity?: { fingerprint?: string, compactSignature?: string } }} handle
 * @returns {unknown}
 */
const artifactRetentionKey = (handle) => {
  const identity = handle?.identity;
  return identity?.fingerprint
    ? `${identity.fingerprint}\u0000${identity.compactSignature ?? ''}`
    : handle;
};

const checkpointSvgBytes = (checkpoint) => (
  Array.isArray(checkpoint?.results)
    ? checkpoint.results.reduce(
        (total, result) => total + String(result?.content || '').length * 2,
        0
      )
    : 0
);

/**
 * @typedef {Record<string, any>} HistoryIntent
 *   The editable intent (config, files, ui, features, ...); the snapshot service builds and applies it.
 * @typedef {Record<string, any>} HistoryCheckpoint
 *   An artifact checkpoint: intent plus the generated Result domains.
 * @typedef {{
 *   label: string,
 *   before: any,
 *   closed: boolean,
 *   source: string,
 *   owner?: any,
 *   deferAdapterCommit?: boolean,
 *   beforeSignature?: string,
 *   beforeFileIds?: Set<any>,
 *   beforeIntentCheckpoint?: any
 * }} HistoryTransaction
 *   An open Undo transaction: an intent edit (`begin`), an artifact checkpoint, or an artifact replacement.
 * @typedef {{ path: string[], before: any, after: any, beforeHas: boolean, afterHas: boolean }} HistoryIntentChange
 *   One changed leaf of an intent; `path[0]` is the domain.
 * @typedef {{
 *   retainedBytes: number,
 *   fileIds?: readonly string[],
 *   identity?: { fingerprint?: string, compactSignature?: string },
 *   [key: string]: any
 * }} HistoryArtifactHandle
 *   A generated-artifact capture; History reads only these fields. Handles with the same
 *   fingerprint and compact signature name one artifact, whose `retainedBytes` count once.
 * @typedef {{ status: string, reason: string }} HistoryBusy
 *   The Session availability that rejects an edit, Undo, or Redo.
 * @typedef {object} HistoryFileRetention
 *   The two functions of the History file store that History calls (R13).
 * @property {(fileIds: Set<string>) => number} estimateBytes
 *   Bytes of the retained files named by `fileIds`.
 * @property {(fileIds: Set<string>) => void} retainOnly
 *   Drops every stored file that `fileIds` does not name.
 *
 * @typedef {object} HistoryManagerOptions
 * @property {() => HistoryIntent | Promise<HistoryIntent>} buildIntent
 *   The snapshot service's capture of the current intent.
 * @property {(intent: HistoryIntent, context: { changes: HistoryIntentChange[], direction: 'undo' | 'redo' }) => unknown} applyIntent
 *   The snapshot service's restore of an intent; History awaits it.
 * @property {() => HistoryCheckpoint | Promise<HistoryCheckpoint>} buildCheckpoint
 *   The snapshot service's capture of an artifact checkpoint.
 * @property {(checkpoint: HistoryCheckpoint) => unknown} applyCheckpoint
 *   The snapshot service's restore of an artifact checkpoint; History awaits it.
 * @property {(() => HistoryArtifactHandle | Promise<HistoryArtifactHandle>) | null} [captureGeneratedArtifactHandle]
 *   The snapshot service's capture of the generated artifact for a replacement step.
 * @property {((handle: HistoryArtifactHandle, options: { clearFailedGeneratePresentation: boolean }) => unknown) | null} [restoreGeneratedArtifactHandle]
 *   The snapshot service's restore of a captured handle; History awaits it.
 * @property {((before: HistoryArtifactHandle, after: HistoryArtifactHandle) => boolean) | null} [compareGeneratedArtifactHandles]
 *   Whether two handles name the same artifact.
 * @property {(value: unknown) => string} [signatureFor]
 *   The signature that detects a change; JSON by default.
 * @property {HistoryFileRetention | null} [fileStore]
 *   The file store that keeps the files that entries name (F-05).
 * @property {(() => Iterable<string>) | null} [collectCurrentFileIds]
 *   The ids of the files bound now, from the snapshot service.
 * @property {number} [maxActions] Undo entries kept.
 * @property {number} [maxBytes] Retained bytes allowed before the oldest entries drop.
 * @property {<T>(value: T) => { value: T }} [makeRef]
 *   Creates a reactive box (Vue `ref` in the app); a plain box by default.
 * @property {(operation?: string) => HistoryBusy | null} [mutationAvailability]
 *   The Session lifecycle's busy answer (`'history'` words Undo and Redo, D-28).
 * @property {(step: { type: string, changes?: unknown }) => HistoryBusy | null} [stepAvailability]
 *   The composition root's busy answer for one Undo or Redo step (E1: a step
 *   that switches the diagram mode waits as the mode buttons do).
 */

/** @param {HistoryManagerOptions} options */
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
  makeRef = makeBox,
  mutationAvailability = () => null,
  stepAvailability = () => null
} = /** @type {HistoryManagerOptions} */ ({})) => {
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
  // Transaction open/close notifies mutationPending() without a document revision.
  const transactionRevision = makeRef(0);
  const restoring = makeRef(false);
  // Undo and Redo calls from their request to their end, the settle of an
  // open intent before `restoring` included (E1: the mode switch waits for them).
  const traversals = makeRef(0);
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
  // OV-112: adjacent artifact replacements name one artifact (entry k's `after`
  // is entry k+1's `before`), so its bytes count once while any live entry
  // names it. Keyed by `artifactRetentionKey`.
  /** @type {Map<unknown, { bytes: number, references: number }>} */
  const retainedArtifacts = new Map();
  /** @type {HistoryTransaction | null} */
  let activeTransaction = null;
  /** @type {HistoryTransaction | null} */
  let activeCheckpoint = null;
  /** @type {HistoryTransaction | null} */
  let activeReplacement = null;

  const touch = () => {
    revision.value += 1;
  };
  const touchTransaction = () => {
    transactionRevision.value += 1;
  };
  const intentCommitOwnedByAction = () => Boolean(
    activeTransaction && !activeTransaction.closed && activeTransaction.deferAdapterCommit
  );

  // D-28 (PD-OI-082): an open artifact replacement or checkpoint, or an action
  // that still owns its intent commit, owns the current state. Undo, Redo, and
  // their buttons and shortcuts share this one reactive availability; the
  // Session owner words a Generate or diagram update in progress.
  const historyAvailability = () => {
    void transactionRevision.value;
    const busy = mutationAvailability();
    if (busy) return busy;
    if (!activeReplacement && !activeCheckpoint && !intentCommitOwnedByAction()) return null;
    return mutationAvailability('history') || HISTORY_EDIT_BUSY;
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

  /**
   * The handles whose bytes an entry names in addition to its own `byteSize`.
   * @param {any} entry
   * @returns {HistoryArtifactHandle[]}
   */
  const entryArtifactHandles = (entry) => (
    entry?.type === 'artifact-replacement' ? [entry.before, entry.after] : []
  );

  /** @param {HistoryArtifactHandle} handle */
  const retainArtifactBytes = (handle) => {
    const key = artifactRetentionKey(handle);
    const bytes = Number(handle?.retainedBytes) || 0;
    const retained = retainedArtifacts.get(key);
    if (!retained) {
      retainedArtifacts.set(key, { bytes, references: 1 });
      totalEntryBytes += bytes;
      return;
    }
    retained.references += 1;
    if (bytes > retained.bytes) {
      totalEntryBytes += bytes - retained.bytes;
      retained.bytes = bytes;
    }
  };

  /** @param {HistoryArtifactHandle} handle */
  const releaseArtifactBytes = (handle) => {
    const key = artifactRetentionKey(handle);
    const retained = retainedArtifacts.get(key);
    if (!retained) return;
    retained.references -= 1;
    if (retained.references > 0) return;
    retainedArtifacts.delete(key);
    totalEntryBytes = Math.max(0, totalEntryBytes - retained.bytes);
  };

  const addEntryBytes = (entry) => {
    totalEntryBytes += Number(entry?.byteSize) || 0;
    entryArtifactHandles(entry).forEach(retainArtifactBytes);
  };

  const removeEntryBytes = (entry) => {
    totalEntryBytes = Math.max(0, totalEntryBytes - (Number(entry?.byteSize) || 0));
    entryArtifactHandles(entry).forEach(releaseArtifactBytes);
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

  // UJ-09: the History position of the last Save or Load. Work is unsaved
  // when an Undo or Redo step was added or traversed since then, or when a
  // control's transaction is still open.
  /** @param {any[]} stack */
  const topEntry = (stack) => stack[stack.length - 1] || null;
  let savePoint = { undo: null, redo: null };
  const markSavePoint = () => {
    savePoint = { undo: topEntry(undoStack), redo: topEntry(redoStack) };
  };
  const hasChangesSinceSavePoint = () => Boolean(activeTransaction && !activeTransaction.closed)
    || topEntry(undoStack) !== savePoint.undo
    || topEntry(redoStack) !== savePoint.redo;

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

  // R11: a transaction belongs to one owner (a control or a gesture). An
  // ownerless action, the same owner, or an action that owns the open commit
  // joins it; another owner first settles it.
  const begin = async (label = 'Edit', options = {}) => {
    if (mutationAvailability()) return null;
    if (restoring.value || capturing.value) return null;
    const open = activeTransaction && !activeTransaction.closed ? activeTransaction : null;
    if (open) {
      if (
        options.owner === undefined
        || open.owner === options.owner
        || open.deferAdapterCommit
      ) return open;
      await settlePendingIntent();
      if (activeTransaction && !activeTransaction.closed) return activeTransaction;
    }

    const before = await captureIntent();
    if (mutationAvailability()) return null;
    emitHistoryDiagnostic({ type: 'begin', scope: 'intent', label });
    /** @type {HistoryTransaction} */
    const tx = {
      label,
      before: before.intent,
      beforeSignature: before.signature,
      beforeFileIds: before.fileIds,
      closed: false,
      source: options.source || '',
      owner: options.owner
    };
    activeTransaction = tx;
    touchTransaction();
    return tx;
  };

  const cancel = (transaction) => {
    if (!transaction) return;
    transaction.closed = true;
    if (activeTransaction === transaction) { activeTransaction = null; touchTransaction(); }
  };

  const initializeIntentBaseline = async (_label = 'Intent baseline', { isCurrent = () => true } = {}) => {
    if (restoring.value) return false;
    if (activeTransaction && !activeTransaction.closed) cancel(activeTransaction);
    const intentRecord = await captureIntent();
    if (!isCurrent()) return false;
    clearStack(undoStack);
    clearRedo();
    markSavePoint();
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
    // A concurrent settlement may have committed it during the capture.
    if (transaction.closed) return false;
    const busy = mutationAvailability();
    if (busy) { cancel(transaction); return busy; }
    transaction.closed = true;
    if (activeTransaction === transaction) { activeTransaction = null; touchTransaction(); }
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

  // R11: the one path that commits the open intent before another owner, a
  // checkpoint, an artifact replacement, a command, Undo, or Redo takes the
  // state. The input adapter skips its own commit of a settled transaction.
  const settlePendingIntent = async () => {
    const pending = activeTransaction;
    if (!pending || pending.closed) return;
    const previousDeferred = Boolean(pending.deferAdapterCommit);
    pending.deferAdapterCommit = true;
    try {
      await commit(pending);
    } catch (error) {
      if (!pending.closed) pending.deferAdapterCommit = previousDeferred;
      throw error;
    }
  };

  const runUndoable = async (label, fn, options = {}) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value || activeCheckpoint) return fn();

    const usesActiveTransaction = Boolean(activeTransaction && !activeTransaction.closed);
    const tx = usesActiveTransaction ? activeTransaction : await begin(label, options);
    if (tx) { tx.deferAdapterCommit = true; touchTransaction(); }

    try {
      const busy = mutationAvailability();
      if (busy) { cancel(tx); return busy; }
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
    return settlePendingIntent().then(captureBefore);
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
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value || activeCheckpoint) return fn();
    const execute = async (tx) => {
      const busy = mutationAvailability();
      if (busy) { cancel(tx); return busy; }
      activeCheckpoint = tx;
      touchTransaction();
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
      } finally {
        activeCheckpoint = null;
        touchTransaction();
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
      const handle = captureGeneratedArtifactHandle();
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

  const captureReplacementIntentCheckpoint = async (options, phase) => {
    if (typeof options?.captureIntentCheckpoint !== 'function') {
      return { enabled: false, value: undefined };
    }
    return {
      enabled: true,
      value: await options.captureIntentCheckpoint({ phase })
    };
  };

  const restoreReplacementIntentCheckpoint = async (checkpoint, restore) => {
    if (!checkpoint?.enabled || typeof restore !== 'function') return;
    restoring.value = true;
    try {
      await restore(checkpoint.value);
    } finally {
      restoring.value = false;
    }
  };

  const beginArtifactReplacement = async (label = 'Generate diagram', options = {}) => {
    if (restoring.value || capturing.value) return null;
    const capturesIntent = typeof options.captureIntentCheckpoint === 'function';
    const restoresIntent = typeof options.restoreIntentCheckpoint === 'function';
    if (capturesIntent !== restoresIntent) {
      throw new Error(
        'Artifact replacement intent checkpoints require capture and restore hooks.'
      );
    }
    await settlePendingIntent();
    notifyCheckpointCapture(options, 'before-start');
    recordSessionLifecycleEvent('history.before-capture-started');
    const before = await captureArtifactHandle('before');
    recordSessionLifecycleEvent('history.before-capture-completed');
    const beforeIntentCheckpoint = await captureReplacementIntentCheckpoint(
      options,
      'before'
    );
    notifyCheckpointCapture(options, 'before-end');
    emitHistoryDiagnostic({ type: 'begin', scope: 'artifact-replacement', label });
    return {
      label,
      before,
      beforeIntentCheckpoint,
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
    const afterIntentCheckpoint = await captureReplacementIntentCheckpoint(
      options,
      'after'
    );
    transaction.closed = true;
    currentFileIds = artifactHandleFileIds(after);
    clearCurrentCheckpoint();

    const intentChanged = transaction.beforeIntentCheckpoint.enabled
      && afterIntentCheckpoint.enabled
      && (typeof options.compareIntentCheckpoints === 'function'
        ? !options.compareIntentCheckpoints(
            transaction.beforeIntentCheckpoint.value,
            afterIntentCheckpoint.value
          )
        : !sameJsonValue(
            transaction.beforeIntentCheckpoint.value,
            afterIntentCheckpoint.value
          ));
    if (sameArtifactHandles(transaction.before, after) && !intentChanged) {
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
    if (intentChanged) {
      collectHistoryFileIds(transaction.beforeIntentCheckpoint.value, fileIds);
      collectHistoryFileIds(afterIntentCheckpoint.value, fileIds);
    }
    const intentCheckpointBytes = intentChanged
      ? (JSON.stringify([
          transaction.beforeIntentCheckpoint.value,
          afterIntentCheckpoint.value
        ]) || '').length * 2
      : 0;
    const namedBytes = intentCheckpointBytes
      + (Number(transaction.before.retainedBytes) || 0)
      + (Number(after.retainedBytes) || 0);
    const entry = {
      type: 'artifact-replacement',
      label: transaction.label || options.label || 'Generate diagram',
      before: transaction.before,
      after,
      // The entry's own bytes; `addEntryBytes` counts `before` and `after` once per artifact.
      byteSize: intentCheckpointBytes,
      fileIds,
      ...(intentChanged ? {
        intentCheckpoint: {
          before: transaction.beforeIntentCheckpoint,
          after: afterIntentCheckpoint,
          restore: options.restoreIntentCheckpoint
        }
      } : {})
    };
    diagnostics.byteEstimateComputations += 1;
    diagnostics.artifactReplacementHistoryEntryCount += 1;
    recordStructuralMetric('historyReplacementCount', 1);
    emitHistoryDiagnostic({
      type: 'size',
      scope: 'artifact-replacement',
      bytes: namedBytes
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
      bytes: namedBytes
    });
    notifyCheckpointCapture(options, 'after-end');
    return true;
  };

  const restoreArtifactHandle = async (handle, { refreshIntent = true } = {}) => {
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
    if (refreshIntent) await refreshCurrentIntent();
  };

  const runUndoableArtifactReplacement = async (label, fn, options = {}) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value) return fn(null);
    const transaction = await beginArtifactReplacement(label, options);
    const busy = mutationAvailability();
    if (busy) { cancel(transaction); return busy; }
    activeReplacement = transaction;
    touchTransaction();
    try {
      const result = await fn(transaction?.before || null);
      if (
        typeof options.shouldCommit === 'function'
        && !options.shouldCommit(result)
      ) {
        if (transaction) transaction.closed = true;
        releaseUnreferencedFiles();
        touch();
        return result;
      }
      recordSessionLifecycleEvent('history.finalization-started', { label });
      await commitAppliedArtifactReplacement(transaction, options);
      recordSessionLifecycleEvent('history.finalization-completed', { label });
      return result;
    } catch (error) {
      try {
        if (transaction) {
          transaction.closed = true;
          if (typeof options.restoreAppliedArtifact === 'function') {
            restoring.value = true;
            try {
              await options.restoreAppliedArtifact(transaction.before);
            } finally {
              restoring.value = false;
            }
            clearCurrentCheckpoint();
          } else {
            recordStructuralMetric('generatedArtifactRollbackCount', 1);
            await restoreArtifactHandle(transaction.before, { refreshIntent: false });
          }
          await restoreReplacementIntentCheckpoint(
            transaction.beforeIntentCheckpoint,
            options.restoreIntentCheckpoint
          );
          await refreshCurrentIntent();
        }
        releaseUnreferencedFiles();
      } catch (_) {
        if (error && typeof error === 'object') error.artifactRestoreFailed = true;
      }
      throw error;
    } finally {
      if (activeReplacement === transaction) { activeReplacement = null; touchTransaction(); }
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
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof buildCommand !== 'function') return false;
    if (!restoring.value && !capturing.value) await settlePendingIntent();
    const preparedCommand = await buildCommand();
    const busy = mutationAvailability();
    if (busy) return busy;
    const command = normalizeCommand(label, preparedCommand);
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

  const applyArtifactReplacementEntry = async (entry, direction) => {
    const handle = direction === 'undo' ? entry.before : entry.after;
    if (!entry.intentCheckpoint) {
      await restoreArtifactHandle(handle);
      return;
    }
    await restoreArtifactHandle(handle, { refreshIntent: false });
    const checkpoint = direction === 'undo'
      ? entry.intentCheckpoint.before
      : entry.intentCheckpoint.after;
    await restoreReplacementIntentCheckpoint(
      checkpoint,
      entry.intentCheckpoint.restore
    );
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

  // SE-04: the focused control's open intent is settled first, so an Undo is
  // never recorded afterward as a new edit that clears Redo.
  const settleBeforeTraversal = async () => {
    const busy = historyAvailability();
    if (busy || restoring.value) return busy;
    await settlePendingIntent();
    return historyAvailability();
  };

  /** @param {() => Promise<unknown>} step */
  const traverse = async (step) => {
    traversals.value += 1;
    try {
      return await step();
    } finally {
      traversals.value -= 1;
    }
  };

  const undo = () => traverse(async () => {
    const busy = await settleBeforeTraversal();
    if (busy) return busy;
    if (restoring.value || undoStack.length === 0) return false;
    const entry = undoStack[undoStack.length - 1];
    const stepBusy = stepAvailability(entry);
    if (stepBusy) return stepBusy;
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
      await applyArtifactReplacementEntry(entry, 'undo');
    } else {
      await applyIntentEntry(entry, 'undo');
    }
    undoStack.pop();
    redoStack.push(entry);
    enforceLimits();
    touch();
    return true;
  });

  const redo = () => traverse(async () => {
    const busy = await settleBeforeTraversal();
    if (busy) return busy;
    if (restoring.value || redoStack.length === 0) return false;
    const entry = redoStack[redoStack.length - 1];
    const stepBusy = stepAvailability(entry);
    if (stepBusy) return stepBusy;
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
      await applyArtifactReplacementEntry(entry, 'redo');
    } else {
      await applyIntentEntry(entry, 'redo');
    }
    redoStack.pop();
    undoStack.push(entry);
    enforceLimits();
    touch();
    return true;
  });

  const canUndo = () => undoStack.length > 0 && !restoring.value && !historyAvailability();
  const canRedo = () => redoStack.length > 0 && !restoring.value && !historyAvailability();
  const undoLabel = () => (canUndo() ? undoStack[undoStack.length - 1].label : '');
  const redoLabel = () => (canRedo() ? redoStack[redoStack.length - 1].label : '');

  return {
    mutationPending: () => {
      void transactionRevision.value;
      return Boolean(restoring.value || capturing.value || intentCommitOwnedByAction());
    },
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
    hasChangesSinceSavePoint,
    historyAvailability,
    historyLimitMessage,
    traversalPending: () => traversals.value > 0,
    initializeIntentBaseline,
    markSavePoint,
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
