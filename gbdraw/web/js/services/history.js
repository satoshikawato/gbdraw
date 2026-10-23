// @ts-check
import { collectHistoryFileIds } from './history-files.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

// tests/web/history-inputs.test.mjs loads this module with its own imports
// only, so it declares the key marker of utils/i18n.js itself.
/** @param {string} text @returns {string} */
const msg = (text) => text;

const DEFAULT_MAX_ACTIONS = 30;
const DEFAULT_MAX_BYTES = 200 * 1024 * 1024;
// The Session availability wording for an edit still being applied.
const HISTORY_EDIT_BUSY = Object.freeze({
  status: 'busy',
  reason: msg('Applying an edit. Retry after the edit finishes.')
});

/**
 * @template T
 * @param {T} value
 * @returns {{ value: T }}
 */
const makeBox = (value) => ({ value });
/**
 * @param {Record<string, unknown>} value
 * @param {string} key
 */
const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);
/**
 * @template T
 * @param {T | PromiseLike<T>} value
 * @returns {value is PromiseLike<T>}
 */
const isPromiseLike = (value) => (
  Boolean(value) && typeof (/** @type {{ then?: unknown }} */ (value)).then === 'function'
);

/**
 * A deep copy of plain data, with the same shape.
 * @template T
 * @param {T} value
 * @returns {T}
 */
const cloneIntentValue = (value) => {
  if (Array.isArray(value)) return /** @type {T} */ (value.map((entry) => cloneIntentValue(entry)));
  if (!isPlainObject(value)) return value;
  return /** @type {T} */ (Object.fromEntries(
    Object.entries(value).map(([key, entry]) => [key, cloneIntentValue(entry)])
  ));
};

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const sameJsonValue = (left, right) => {
  if (Object.is(left, right)) return true;
  if (typeof left !== typeof right) return false;
  if (left === null || right === null) return false;
  if (typeof left !== 'object') return false;
  return JSON.stringify(left) === JSON.stringify(right);
};

/**
 * @param {unknown} before
 * @param {unknown} after
 * @param {string[]} [path]
 * @param {HistoryIntentChange[]} [changes]
 * @returns {HistoryIntentChange[]}
 */
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

/**
 * @param {HistoryIntent | null | undefined} source
 * @param {readonly HistoryIntentChange[]} changes
 * @param {'undo' | 'redo'} direction
 * @returns {HistoryIntent}
 */
const applyIntentPatch = (source, changes, direction) => {
  /** @type {HistoryIntent} */
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

/**
 * The value at `path` in an intent, and whether it is present.
 * @param {HistoryIntent | null} intent
 * @param {readonly string[]} path
 * @returns {{ has: boolean, value: unknown }}
 */
const intentValueAt = (intent, path) => {
  /** @type {unknown} */
  let at = intent;
  for (const key of path) {
    if (!isPlainObject(at) || !hasOwn(at, key)) return { has: false, value: undefined };
    at = /** @type {Record<string, unknown>} */ (at)[key];
  }
  return { has: true, value: at };
};

/**
 * The changes of a step that `intent` still holds as the step left them
 * (`direction` leaves that side): what an Undo or Redo restores of a step
 * whose values a later step elsewhere may have changed again (OIC-030).
 * `editedPartOf` gives the part of a value a later edit changes (the intent
 * owner's `editedPartOf` port).
 * @param {HistoryIntent | null} intent
 * @param {readonly HistoryIntentChange[]} changes
 * @param {'undo' | 'redo'} direction
 * @param {IntentEditedPart} [editedPartOf]
 */
const changesAsLeft = (intent, changes, direction, editedPartOf = (_path, value) => value) => {
  const leaving = direction === 'undo' ? 'after' : 'before';
  return changes.filter((change) => {
    const { has, value } = intentValueAt(intent, change.path);
    return has === change[`${leaving}Has`]
      && (!has || sameJsonValue(editedPartOf(change.path, value), editedPartOf(change.path, change[leaving])));
  });
};

/** @param {Record<string, unknown>} event */
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

/** @param {HistoryCheckpoint | null | undefined} checkpoint */
const checkpointSvgBytes = (checkpoint) => (
  Array.isArray(checkpoint?.results)
    ? checkpoint.results.reduce(
        (/** @type {number} */ total, /** @type {{ content?: unknown } | null | undefined} */ result) =>
          total + String(result?.content || '').length * 2,
        0
      )
    : 0
);

/**
 * What one History step captures about the intent: the intent, its signature, and the files it names.
 * @typedef {{ intent: HistoryIntent, signature: string, fileIds: Set<string> }} IntentRecord
 * @typedef {{ checkpoint: HistoryCheckpoint, signature: string, byteSize: number, fileIds: Set<string>, drawn?: boolean }} CheckpointRecord
 *   `drawn`: the step's caller said whether the figure is the one its state
 *   draws, so that restoring it needs no automatic rerender (D-24 (OV-731)).
 * @typedef {{ enabled: boolean, value: unknown }} ReplacementIntentCheckpoint
 *   The intent a replacement step also restores (`value` is the caller's own capture).
 * @typedef {{
 *   type: 'intent', label: string, changes: HistoryIntentChange[], byteSize: number, fileIds: Set<string>, sequence?: number,
 *   session?: IntentEntry
 * }} IntentEntry
 * @typedef {{
 *   type: 'checkpoint', label: string, before: CheckpointRecord, after: CheckpointRecord, byteSize: number,
 *   fileIds: Set<string>, sequence?: number, changes?: HistoryIntentChange[], session?: IntentEntry
 * }} CheckpointEntry
 *   `changes`: the intent values the step changed, kept by a part of a project step (`takeEntry`).
 * @typedef {{
 *   type: 'artifact-replacement',
 *   label: string,
 *   before: HistoryArtifactHandle,
 *   after: HistoryArtifactHandle,
 *   byteSize: number,
 *   fileIds: Set<string>,
 *   sequence?: number,
 *   session?: IntentEntry,
 *   intentCheckpoint?: {
 *     before: ReplacementIntentCheckpoint,
 *     after: ReplacementIntentCheckpoint,
 *     restore: HistoryStepBase['restoreIntentCheckpoint']
 *   }
 * }} ArtifactReplacementEntry
 * @typedef {{
 *   type: 'command', label: string, apply: () => unknown, revert: () => unknown, byteSize: number, fileIds: Set<string>,
 *   sequence?: number, session?: IntentEntry
 * }} CommandEntry
 * @typedef {IntentEntry | CheckpointEntry | ArtifactReplacementEntry | CommandEntry} HistoryEntry
 *   One Undo or Redo step. `sequence`: the order in which the router's steps ended (D-15-52).
 *   `session`: the Session-level values the step changed with this drawing's (D-15-53), restored with it.
 * @typedef {{ noop?: boolean, label?: string, apply: () => unknown, revert: () => unknown, estimateBytes?: () => number, metadata?: unknown }} UndoableCommand
 *   A reversible command; `apply` or `revert` returning `false` means it did not run.
 * @typedef {{
 *   capture(context: { phase: 'before' | 'after' }): unknown,
 *   restore(checkpoint: unknown): unknown
 * }} IntentCheckpointHooks
 *   The two intent checkpoint hooks a caller may also pass as `null`.
 * @typedef {{
 *   onCheckpointCapture?: ((event: { phase: string, diagnostics: Record<string, number> }) => void) | null,
 *   label?: string,
 *   source?: string,
 *   owner?: unknown,
 *   joinsOpenStep?: boolean,
 *   captureIntentCheckpoint?: IntentCheckpointHooks['capture'] | null,
 *   restoreIntentCheckpoint?: IntentCheckpointHooks['restore'] | null,
 *   compareIntentCheckpoints?(before: unknown, after: unknown): boolean,
 *   restoreAppliedArtifact?: (before: HistoryArtifactHandle) => unknown,
 *   drawn?: { before: () => boolean, change: boolean },
 *   takeEntry?: (entry: HistoryEntry | null) => void
 * }} HistoryStepBase
 *   What a step reads from its caller; a hook that is not a function is not called. The intent
 *   checkpoint is the caller's own value, so the three hooks take and return it opaque (they are
 *   declared as methods); a caller names its type on its own parameter. A checkpoint step also
 *   reads `drawn`, whether the figure of each side is the one its state draws (`CheckpointRecord`):
 *   `before()` as the before side is captured, and `change`, whether the step's own change shows
 *   as its state draws it, so that the after side is drawn when the before side is. `takeEntry`:
 *   the router's project step takes the step's entry (or `null`) in place of this History's Undo stack.
 */
/**
 * @template [R=unknown]
 * @typedef {HistoryStepBase & { shouldCommit?: (result: R) => boolean }} HistoryStepOptions
 *   A step that runs `fn` commits only when `shouldCommit(result)` allows it.
 */

/**
 * @typedef {Record<string, any>} HistoryIntent
 *   The editable intent (config, files, ui, features, ...); the snapshot service builds and applies it.
 * @typedef {Record<string, any>} HistoryCheckpoint
 *   An artifact checkpoint: intent plus the generated Result domains.
 * @typedef {{ intent: HistoryIntent, signature: string }} SessionMark
 *   The router's record of the Session-level values as a step opened (D-15-53).
 * @typedef {{
 *   label: string,
 *   before: any,
 *   closed: boolean,
 *   committedAt?: number,
 *   actions?: { label: string, stepLabel?: string }[],
 *   source: string,
 *   owner?: any,
 *   deferAdapterCommit?: boolean,
 *   serial?: number,
 *   beforeSignature?: string,
 *   beforeFileIds?: Set<any>,
 *   beforeIntentCheckpoint?: any,
 *   beforeIntent?: HistoryIntent | null,
 *   takeEntry?: HistoryStepBase['takeEntry'],
 *   sessionMark?: SessionMark | null
 * }} HistoryTransaction
 *   An open Undo transaction: an intent edit (`begin`), an artifact checkpoint, or an artifact replacement.
 *   `committedAt`: the intent generation at which a commit closed it. `actions`: the `runUndoable`
 *   calls that ran on it, each with the step name it started with. `serial`: an intent step's
 *   number, in the order the steps began (`intentStepSerial`). `sessionMark`: the router's Session
 *   values as the step opened (`stepOpened`).
 * @typedef {{ path: string[], before: any, after: any, beforeHas: boolean, afterHas: boolean }} HistoryIntentChange
 *   One changed leaf of an intent; `path[0]` is the domain.
 * @typedef {(path: readonly string[], value: unknown) => unknown} IntentEditedPart
 *   The part of the intent value at `path` that an edit changes.
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
 * @property {(checkpoint: HistoryCheckpoint, restored: { drawn: boolean }) => unknown} applyCheckpoint
 *   The snapshot service's restore of an artifact checkpoint; History awaits it. `restored.drawn`:
 *   whether its figure is the one its state draws.
 * @property {(() => HistoryArtifactHandle | Promise<HistoryArtifactHandle>) | null} [captureGeneratedArtifactHandle]
 *   The snapshot service's capture of the generated artifact for a replacement step.
 * @property {((handle: HistoryArtifactHandle, options: { clearFailedGeneratePresentation: boolean }) => unknown) | null} [restoreGeneratedArtifactHandle]
 *   The snapshot service's restore of a captured handle; History awaits it.
 * @property {((before: HistoryArtifactHandle, after: HistoryArtifactHandle) => boolean) | null} [compareGeneratedArtifactHandles]
 *   Whether two handles name the same artifact.
 * @property {(value: unknown) => string} [signatureFor]
 *   The signature that detects a change; JSON by default.
 * @property {HistoryFileRetention | null} [fileStore]
 *   The file store that keeps the files that entries name (F-05); read only
 *   when no `budget` is given.
 * @property {(() => Iterable<string>) | null} [collectCurrentFileIds]
 *   The ids of the files bound now, from the snapshot service.
 * @property {number} [maxActions] Undo entries kept.
 * @property {number} [maxBytes] Retained bytes allowed before the oldest entries drop; read
 *   only when no `budget` is given.
 * @property {HistoryBudget} [budget]
 *   The byte budget and file retention this manager shares with the other
 *   drawings' managers (`createHistoryRouter`); its own by default.
 * @property {(message: string) => void} [onHistoryLimit]
 *   Told when the byte limit dropped steps, with the sentence that says so
 *   (English, marked with `msg`); the step-count window (`maxActions`) drops
 *   its oldest step silently.
 * @property {<T>(value: T) => { value: T }} [makeRef]
 *   Creates a reactive box (Vue `ref` in the app); a plain box by default.
 * @property {(operation?: string) => HistoryBusy | null} [mutationAvailability]
 *   The Session lifecycle's busy answer (`'history'` words Undo and Redo, D-28).
 * @property {(step: { type: string, changes?: unknown }) => HistoryBusy | null} [stepAvailability]
 *   The composition root's busy answer for one Undo or Redo step. Undo, Redo,
 *   and their enabled state read it for the top step, so it must read reactive state.
 * @property {readonly string[] | null} [identityMapPath]
 *   The intent path of a map keyed by Result identity (the composition
 *   offsets); `heldIdentities` lists the keys History holds there (OV-110).
 * @property {IntentEditedPart} [editedPartOf]
 *   The part of an intent value that an edit changes, without what the drawn
 *   figure lays out; a project part restores a value whose edited part is
 *   still as its step left it (OV-1655). The whole value by default.
 * @property {() => SessionMark | null} [stepOpened]
 *   The router's Session values as a step opens here; the step carries them to its end.
 * @property {(entry: HistoryEntry | null, label: string, mark: SessionMark | null) => IntentEntry | null} [stepEnded]
 *   The router hears each step that ends on this History's Undo stack, with
 *   its entry or `null` when it changed nothing here (D-15-51, R15-23), and
 *   the Session values it opened with; it answers the Session-level values
 *   the entry's step also changed (D-15-53).
 * @property {(part: IntentEntry, direction: 'undo' | 'redo') => unknown} [applySessionPart]
 *   Restores an entry's `session` part, where still as its step left it.
 */

/**
 * What the shared budget reads from and drops in one History manager.
 * @typedef {object} HistoryBudgetMember
 * @property {() => number} entryBytes The bytes its entries retain.
 * @property {() => Set<string>} retainedFileIds The files its current intent and entries name.
 * @property {(kept: object | null) => 'undo' | 'redo' | null} dropOldestStep
 *   Drops its oldest Undo entry, or its Redo entries when no Undo entry is left,
 *   and says which; `null` when it holds none. It keeps `kept`, the step just
 *   recorded, undone, or redone (OV-1260).
 * @property {(message: string) => void} noteEvicted
 *   Tells that the byte limit dropped its steps, with the sentence that says so.
 *
 * The byte budget and file retention of every drawing's History (OV-112).
 * @typedef {object} HistoryBudget
 * @property {(member: HistoryBudgetMember) => void} join
 * @property {(kept?: object | null) => void} enforce
 *   Drops steps until the retained bytes fit, oldest first, from the members in
 *   `evictionOrder`, keeping `kept` (OV-1260); then releases the files no member names.
 * @property {() => void} releaseUnreferencedFiles
 */

/**
 * @typedef {object} HistoryBudgetOptions
 * @property {number} [maxBytes] Retained bytes allowed across the members.
 * @property {(members: readonly HistoryBudgetMember[]) => HistoryBudgetMember[]} [evictionOrder]
 *   The order in which the members lose steps; the order they joined by default.
 * @property {HistoryFileRetention['estimateBytes'] | null} [estimateFileBytes] The file store's estimate.
 * @property {HistoryFileRetention['retainOnly'] | null} [retainOnlyFiles] The file store's release.
 */
/**
 * @param {HistoryBudgetOptions} [options]
 * @returns {HistoryBudget}
 */
export const createHistoryBudget = ({
  maxBytes = DEFAULT_MAX_BYTES, estimateFileBytes = null, retainOnlyFiles = null, evictionOrder = (members) => [...members]
} = {}) => {
  /** @type {Set<HistoryBudgetMember>} */
  const members = new Set();
  const retainedFileIds = () => {
    /** @type {Set<string>} */
    const ids = new Set();
    members.forEach((member) => member.retainedFileIds().forEach((id) => ids.add(id)));
    return ids;
  };
  const overBudget = () => {
    let bytes = 0;
    members.forEach((member) => { bytes += member.entryBytes(); });
    if (estimateFileBytes) bytes += estimateFileBytes(retainedFileIds());
    return bytes > maxBytes;
  };
  const releaseUnreferencedFiles = () => {
    if (retainOnlyFiles) retainOnlyFiles(retainedFileIds());
  };
  // OV-1260, OV-1365: every byte-limit drop is announced through the one
  // notice path, in English marked with `msg`: the notice translates it where
  // it shows it.
  /** @param {object | null} [kept] */
  const enforce = (kept = null) => {
    evictionOrder([...members]).forEach((candidate) => {
      let undoDropped = false;
      let redoDropped = false;
      while (overBudget()) {
        const dropped = candidate.dropOldestStep(kept);
        if (!dropped) break;
        if (dropped === 'undo') undoDropped = true;
        else redoDropped = true;
      }
      const message = undoDropped && redoDropped ? msg('Older undo history and the redo history were discarded to stay within the history limit.')
        : undoDropped ? msg('Older undo history was discarded to stay within the history limit.')
          : redoDropped ? msg('Redo history was discarded to stay within the history limit.') : '';
      if (message) candidate.noteEvicted(message);
    });
    releaseUnreferencedFiles();
  };
  return { join: (member) => { members.add(member); }, enforce, releaseUnreferencedFiles };
};

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
  budget = createHistoryBudget({
    maxBytes, estimateFileBytes: fileStore?.estimateBytes, retainOnlyFiles: fileStore?.retainOnly
  }),
  onHistoryLimit = () => {},
  makeRef = makeBox,
  mutationAvailability = () => null,
  stepAvailability = () => null,
  identityMapPath = null,
  editedPartOf = (_path, value) => value,
  stepOpened = () => null,
  stepEnded = () => null,
  applySessionPart = () => {}
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

  /** @type {HistoryEntry[]} */
  const undoStack = [];
  /** @type {HistoryEntry[]} */
  const redoStack = [];
  const revision = makeRef(0);
  // Transaction open/close notifies mutationPending() without a document revision.
  const transactionRevision = makeRef(0);
  const restoring = makeRef(false);
  // Undo and Redo calls from their request to their end, the settle of an
  // open intent before `restoring` included (E1: the mode switch waits for them).
  const traversals = makeRef(0);
  const capturing = makeRef(false);
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
  /** @type {HistoryIntent | null} */
  let currentIntent = null;
  let currentIntentSignature = '';
  /** @type {HistoryCheckpoint | null} */
  let currentCheckpoint = null;
  let currentCheckpointSignature = '';
  /** @type {Set<string>} */
  let currentFileIds = new Set();
  // OV-463: advanced where the intent changes outside the intent steps: an
  // artifact replacement commits or ends without committing (the last
  // recorded intent is not refreshed) or a baseline replaces the stacks. An
  // action's later write is recorded from the last recorded intent only
  // within the generation of its step.
  let intentGeneration = 0;
  // D-40: the number of the last intent step begun; a redraw requested at
  // serial n is drawn for the steps up to n (`adoptDrawnWrite`). A step
  // reserves its number when `openStep` starts (some numbers stay unused).
  let intentStepSerial = 0;
  let totalEntryBytes = 0;
  /** @type {(changes: readonly HistoryIntentChange[] | null) => void} */
  let settingsStepListener = () => {};
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
  const openTransaction = () => (
    activeTransaction && !activeTransaction.closed ? activeTransaction : null
  );
  const intentCommitOwnedByAction = () => Boolean(openTransaction()?.deferAdapterCommit);

  // OV-110: the keys under `identityMapPath` that the current intent, the open
  // step, and the intent steps of Undo and Redo hold.
  /** @param {readonly HistoryEntry[]} [projectParts] This drawing's parts of the router's project steps. */
  const heldIdentities = (projectParts = []) => {
    /** @type {Set<string>} */
    const ids = new Set();
    if (!identityMapPath) return ids;
    const mapPath = identityMapPath;
    /** @param {unknown} value @param {readonly (string | number)[]} rest */
    const addKeysBelow = (value, rest) => {
      const map = rest.reduce((/** @type {unknown} */ at, key) => (isPlainObject(at) ? /** @type {Record<string, unknown>} */ (at)[key] : undefined), value);
      if (isPlainObject(map)) Object.keys(/** @type {object} */ (map)).forEach((key) => ids.add(key));
    };
    addKeysBelow(currentIntent, mapPath);
    const open = openTransaction();
    if (open) addKeysBelow(open.before, mapPath);
    [...undoStack, ...redoStack, ...projectParts].forEach((entry) => {
      const changes = entry.type === 'intent' || entry.type === 'checkpoint' ? entry.changes : null;
      (changes || []).forEach(({ path, before, after }) => {
        const shared = Math.min(path.length, mapPath.length);
        if (!mapPath.slice(0, shared).every((key, index) => path[index] === key)) return;
        if (path.length > mapPath.length) {
          ids.add(String(path[mapPath.length]));
          return;
        }
        addKeysBelow(before, mapPath.slice(path.length));
        addKeysBelow(after, mapPath.slice(path.length));
      });
    });
    return ids;
  };

  // D-28 (PD-OI-082): an open artifact replacement or checkpoint, or an action
  // that still owns its intent commit, owns the current state. Undo, Redo, and
  // their buttons and shortcuts share this one reactive availability; Undo,
  // Redo, and their enabled state also read the top step's check
  // (`topStepAvailability`, OV-1339). The Session owner words a Generate or
  // diagram update in progress.
  const historyAvailability = () => {
    void transactionRevision.value;
    const busy = mutationAvailability();
    if (busy) return busy;
    if (!activeReplacement && !activeCheckpoint && !intentCommitOwnedByAction()) return null;
    return mutationAvailability('history') || HISTORY_EDIT_BUSY;
  };

  /**
   * @param {unknown} value
   * @param {string} scope
   */
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

  /**
   * @param {string} signature
   * @param {string} scope
   */
  const estimateSignedBytes = (signature, scope) => {
    diagnostics.byteEstimateComputations += 1;
    const bytes = String(signature || '').length * 2;
    emitHistoryDiagnostic({ type: 'size', scope, bytes });
    return bytes;
  };

  /**
   * @param {HistoryStepBase | undefined} options
   * @param {string} phase
   */
  const notifyCheckpointCapture = (options, phase) => {
    const callback = options?.onCheckpointCapture;
    if (typeof callback !== 'function') return;
    try {
      callback({ phase, diagnostics: { ...diagnostics } });
    } catch (_error) {
      // Diagnostic observers cannot own History behavior.
    }
  };

  /** @param {HistoryIntent} intent */
  const collectIntentFileIds = (intent) => {
    const fileIds = collectHistoryFileIds(intent, new Set());
    if (typeof collectCurrentFileIds !== 'function') return fileIds;
    const currentIds = collectCurrentFileIds();
    if (currentIds && typeof currentIds[Symbol.iterator] === 'function') {
      for (const fileId of currentIds) fileIds.add(fileId);
    }
    return fileIds;
  };

  /**
   * @param {HistoryIntent} intent
   * @returns {IntentRecord}
   */
  const finishIntentCapture = (intent) => {
    const signature = computeSignature(intent, 'intent');
    emitHistoryDiagnostic({ type: 'capture', scope: 'intent' });
    return {
      intent,
      signature,
      fileIds: collectIntentFileIds(intent)
    };
  };

  // OV-922: a synchronous capture writes no reactive state; nothing else
  // runs while it does. `capturing` marks a capture only across its await,
  // so a tracked control renders the app for its step, not for each capture.
  /** @returns {IntentRecord | Promise<IntentRecord>} */
  const captureIntent = () => {
    diagnostics.intentBuilds += 1;
    const intent = buildIntent();
    if (!isPromiseLike(intent)) return finishIntentCapture(intent);
    capturing.value = true;
    return Promise.resolve(intent)
      .then(finishIntentCapture)
      .finally(() => {
        capturing.value = false;
      });
  };

  /**
   * @param {HistoryCheckpoint} checkpoint
   * @returns {CheckpointRecord}
   */
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

  /** @param {{ intent: HistoryIntent | null, signature: string, fileIds?: Iterable<string> | null }} record */
  const setCurrentIntent = ({ intent, signature, fileIds }) => {
    currentIntent = intent;
    currentIntentSignature = signature;
    currentFileIds = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
  };

  /** @param {{ checkpoint: HistoryCheckpoint, signature: string }} record */
  const setCurrentCheckpoint = ({ checkpoint, signature }) => {
    currentCheckpoint = checkpoint;
    currentCheckpointSignature = signature;
  };

  const clearCurrentCheckpoint = () => {
    currentCheckpoint = null;
    currentCheckpointSignature = '';
  };

  /**
   * @param {HistoryEntry | null | undefined} entry
   * @returns {Set<string>}
   */
  const entryFileIds = (entry) => (
    entry?.fileIds instanceof Set ? entry.fileIds : new Set()
  );

  const collectRetainedFileIds = () => {
    const ids = new Set(currentFileIds);
    undoStack.forEach((entry) => entryFileIds(entry).forEach((id) => ids.add(id)));
    redoStack.forEach((entry) => entryFileIds(entry).forEach((id) => ids.add(id)));
    return ids;
  };

  const { releaseUnreferencedFiles } = budget;

  /**
   * The handles whose bytes an entry names in addition to its own `byteSize`.
   * @param {HistoryEntry | null | undefined} entry
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

  /** @param {HistoryEntry} entry */
  const addEntryBytes = (entry) => {
    totalEntryBytes += (Number(entry?.byteSize) || 0) + (entry.session?.byteSize ?? 0);
    entryArtifactHandles(entry).forEach(retainArtifactBytes);
  };

  /** @param {HistoryEntry | undefined} entry */
  const removeEntryBytes = (entry) => {
    totalEntryBytes = Math.max(0, totalEntryBytes - (Number(entry?.byteSize) || 0) - (entry?.session?.byteSize ?? 0));
    entryArtifactHandles(entry).forEach(releaseArtifactBytes);
  };

  /** @param {HistoryEntry[]} stack */
  const clearStack = (stack) => {
    stack.forEach(removeEntryBytes);
    stack.splice(0, stack.length);
  };

  /** @param {HistoryEntry} entry */
  const pushUndoEntry = (entry) => {
    undoStack.push(entry);
    addEntryBytes(entry);
  };

  const clearRedo = () => clearStack(redoStack);

  // D-15-51, R15-23: a step's entry goes on this History's Undo stack, with
  // the Session-level values its step also changed (D-15-53), or to the
  // project step that took it (`takeEntry`, Reset all drawings). Either way a
  // new entry clears Redo.
  /**
   * @param {HistoryEntry | null} entry
   * @param {{ label: string, takeEntry?: HistoryStepBase['takeEntry'], sessionMark?: SessionMark | null }} step
   */
  const endStep = (entry, step) => {
    if (step.takeEntry) {
      step.takeEntry(entry);
    } else {
      const session = stepEnded(entry, step.label, step.sessionMark ?? null);
      if (entry && session) entry.session = session;
    }
    if (!entry) return;
    if (!step.takeEntry) pushUndoEntry(entry);
    clearRedo();
  };

  // UJ-09: the History position of the last Save or Load. Work is unsaved
  // when an Undo or Redo step was added or traversed since then, when an
  // open control transaction already changed the intent it began with, or
  // when a value outside History changed (D-07 (b): the Rich Feature Popup
  // preference, which records no step).
  /** @param {HistoryEntry[]} stack */
  const topEntry = (stack) => stack[stack.length - 1] || null;
  /** @type {{ undo: HistoryEntry | null, redo: HistoryEntry | null }} */
  let savePoint = { undo: null, redo: null };
  let changedOutsideHistory = false;
  const markSavePoint = () => {
    savePoint = { undo: topEntry(undoStack), redo: topEntry(redoStack) };
    changedOutsideHistory = false;
  };
  // Records no step and leaves Redo and an open transaction as they are.
  const markChangeOutsideHistory = () => { changedOutsideHistory = true; };
  const hasChangesSinceSavePoint = async () => {
    if (changedOutsideHistory) return true;
    if (topEntry(undoStack) !== savePoint.undo || topEntry(redoStack) !== savePoint.redo) return true;
    const open = openTransaction();
    if (!open) return false;
    // A text control begins on focus; focus alone changes nothing.
    return (await captureIntent()).signature !== open.beforeSignature;
  };

  // OV-1169: the newest step the limit or the budget dropped from Undo.
  let droppedThrough = -1;
  const dropOldestUndo = () => {
    const entry = undoStack.shift();
    droppedThrough = Math.max(droppedThrough, entry?.sequence ?? -1);
    removeEntryBytes(entry);
  };
  /** @type {HistoryBudgetMember} */
  const budgetMember = {
    entryBytes: () => totalEntryBytes,
    retainedFileIds: collectRetainedFileIds,
    // OV-1260: the step just recorded, undone, or redone (`kept`) stays even
    // when it alone is over the byte limit; the oldest other steps drop first.
    dropOldestStep: (kept) => {
      if (undoStack.length > 0 && undoStack[0] !== kept) {
        dropOldestUndo();
        return 'undo';
      }
      const keptRedo = redoStack.filter((entry) => entry === kept);
      if (keptRedo.length === redoStack.length) return null;
      redoStack.filter((entry) => entry !== kept).forEach(removeEntryBytes);
      redoStack.splice(0, redoStack.length, ...keptRedo);
      return 'redo';
    },
    noteEvicted: onHistoryLimit
  };
  budget.join(budgetMember);

  // OV-1365: the step-count window drops silently, as an undo window does;
  // only drops caused by the byte limit are announced.
  /** @param {HistoryEntry | null} [kept] */
  const enforceLimits = (kept = null) => {
    while (undoStack.length > maxActions) dropOldestUndo();
    budget.enforce(kept);
  };

  const captureBaseline = async (_label = 'Baseline') => {
    if (restoring.value) return;
    const checkpointRecord = await captureCheckpoint();
    const intentRecord = await captureIntent();
    clearStack(undoStack);
    clearRedo();
    intentGeneration += 1;
    setCurrentCheckpoint(checkpointRecord);
    setCurrentIntent(intentRecord);
    enforceLimits();
    touch();
  };

  // R11: a transaction belongs to one owner (a control or a gesture). An
  // ownerless action, the same owner, or an action that owns the open commit
  // joins it; another owner first settles it.
  /**
   * @param {string} label
   * @param {HistoryStepBase} options
   * @param {boolean} ownedByAction `runUndoable` opens the step and commits it.
   * @returns {Promise<HistoryTransaction | null>}
   */
  const openStep = async (label, options, ownedByAction) => {
    if (mutationAvailability()) return null;
    if (restoring.value || capturing.value) return null;
    // D-40 (N-A of #1091): numbered as it begins, before any wait below, so a
    // rerender its own edit requests in this task counts it as the requester.
    // A step that joins the open one or opens none leaves an unused number.
    const serial = ++intentStepSerial;
    const open = openTransaction();
    if (open) {
      if (
        options.owner === undefined
        || open.owner === options.owner
        || open.deferAdapterCommit
      ) return open;
      await settlePendingIntent();
      const reopened = openTransaction();
      if (reopened) return reopened;
    }

    // OV-946: a synchronous capture opens the step in this turn, so a
    // control that begins meanwhile finds it and joins it; an asynchronous
    // one holds `capturing` until it ends.
    const capture = captureIntent();
    const before = isPromiseLike(capture) ? await capture : capture;
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
      owner: options.owner,
      deferAdapterCommit: ownedByAction,
      serial,
      takeEntry: options.takeEntry,
      sessionMark: stepOpened()
    };
    activeTransaction = tx;
    touchTransaction();
    return tx;
  };

  /**
   * @param {string} [label]
   * @param {HistoryStepBase} [options]
   */
  const begin = (label = 'Edit', options = {}) => openStep(label, options, false);

  /** @param {HistoryTransaction | null | undefined} transaction */
  const cancel = (transaction) => {
    if (!transaction) return;
    transaction.closed = true;
    if (activeTransaction === transaction) { activeTransaction = null; touchTransaction(); }
  };

  /** @param {string} [_label] @param {{ isCurrent?: () => boolean }} [options] */
  const initializeIntentBaseline = async (_label = 'Intent baseline', { isCurrent = () => true } = {}) => {
    if (restoring.value) return false;
    cancel(openTransaction());
    const intentRecord = await captureIntent();
    if (!isCurrent()) return false;
    clearStack(undoStack);
    clearRedo();
    intentGeneration += 1;
    markSavePoint();
    currentCheckpoint = null;
    currentCheckpointSignature = '';
    setCurrentIntent(intentRecord);
    enforceLimits();
    touch();
    return true;
  };

  /**
   * @param {Pick<HistoryTransaction, 'label' | 'before' | 'beforeSignature' | 'beforeFileIds'>} transaction
   * @param {IntentRecord} afterRecord
   * @param {{ label?: string }} [options]
   * @returns {IntentEntry | null}
   */
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

  /**
   * The bookkeeping of every intent step.
   * @param {Pick<HistoryTransaction, 'label' | 'before' | 'beforeSignature' | 'beforeFileIds' | 'takeEntry' | 'sessionMark'>} from
   *   The transaction, or the last recorded intent when an action records its remainder.
   * @param {Awaited<ReturnType<typeof captureIntent>>} afterRecord
   * @param {{ label?: string }} options
   */
  const recordIntentStep = (from, afterRecord, options) => {
    setCurrentIntent(afterRecord);

    const entry = buildPatchEntry(from, afterRecord, options);
    endStep(entry, { label: entry?.label ?? from.label, takeEntry: from.takeEntry, sessionMark: from.sessionMark });
    if (!entry) {
      emitHistoryDiagnostic({ type: 'commit', scope: 'intent', label: from.label, created: false });
      releaseUnreferencedFiles();
      touch();
      return false;
    }

    emitHistoryDiagnostic({
      type: 'commit',
      scope: 'intent',
      label: entry.label,
      created: true,
      changeCount: entry.changes.length,
      bytes: entry.byteSize
    });
    enforceLimits(entry);
    touch();
    settingsStepListener(entry.changes);
    return true;
  };

  /**
   * @param {HistoryTransaction | null} transaction
   * @param {{ label?: string }} [options]
   * @param {{ label: string, stepLabel?: string } | null} [action]
   *   The `runUndoable` action whose `fn` has ended: its own label, and the name of the step it commits.
   */
  const commitIntent = async (transaction, options = {}, action = null) => {
    if (!transaction) return false;
    // OV-460: a joined action, a drag end, or a settlement may commit an
    // action's transaction while its `fn` awaits. What `fn` wrote afterwards
    // is the action's remainder, recorded as a step of its own. A cancelled
    // transaction has none, an open checkpoint captures it, and a restore
    // owns the stacks while it runs. OV-463: nor is there one under an open
    // artifact replacement or in a later intent generation, where the last
    // recorded intent also lacks what the replacement or the baseline changed.
    const recordsRemainder = () => Boolean(
      action && transaction.committedAt === intentGeneration
      && !activeCheckpoint && !activeReplacement && !restoring.value
    );
    const settled = () => transaction.closed && !recordsRemainder();
    if (settled()) return false;
    const afterRecord = await captureIntent();
    // A concurrent settlement may have committed it during the capture.
    if (settled()) return false;
    const busy = mutationAvailability();
    if (busy) { cancel(transaction); return busy; }
    if (!transaction.closed) {
      transaction.closed = true;
      transaction.committedAt = intentGeneration;
      if (activeTransaction === transaction) { activeTransaction = null; touchTransaction(); }
      return recordIntentStep(
        action ? { ...transaction, label: action.stepLabel ?? transaction.label } : transaction,
        afterRecord,
        options
      );
    }
    if (!action || afterRecord.signature === currentIntentSignature) return false;
    const remainder = {
      label: action.label,
      before: currentIntent,
      beforeSignature: currentIntentSignature,
      beforeFileIds: currentFileIds
    };
    // A transaction opened meanwhile (a control the user is typing in)
    // continues from here, so its step does not repeat the remainder.
    const open = openTransaction();
    if (open) {
      open.before = afterRecord.intent;
      open.beforeSignature = afterRecord.signature;
      open.beforeFileIds = afterRecord.fileIds;
    }
    return recordIntentStep(remainder, afterRecord, options);
  };

  /**
   * @param {HistoryTransaction | null} transaction
   * @param {{ label?: string }} [options]
   */
  const commit = (transaction, options = {}) => commitIntent(transaction, options);

  // R11: the one path that commits the open intent before another owner, a
  // checkpoint, an artifact replacement, a command, Undo, or Redo takes the
  // state. The input adapter skips its own commit of a settled transaction.
  const settlePendingIntent = async () => {
    const pending = openTransaction();
    if (!pending) return;
    const previousDeferred = Boolean(pending.deferAdapterCommit);
    pending.deferAdapterCommit = true;
    try {
      await commit(pending);
    } catch (error) {
      if (!pending.closed) pending.deferAdapterCommit = previousDeferred;
      throw error;
    }
  };

  // The `runUndoable` calls whose `fn` is being called (not awaited) now.
  let actionCallDepth = 0;
  /** @param {() => unknown} fn */
  const callAction = (fn) => {
    actionCallDepth += 1;
    try { return fn(); } finally { actionCallDepth -= 1; }
  };

  /**
   * One undoable step around `fn`, committed when `fn` ends.
   * @param {string} label
   * @param {() => unknown} fn
   * @param {{ joinsOpenStep?: boolean, owner?: unknown, source?: string, label?: string }} [options]
   *   `joinsOpenStep`: the call is part of the action that runs on the open step.
   */
  const runUndoable = async (label, fn, options = {}) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value || activeCheckpoint) return fn();

    const open = openTransaction();
    const tx = open || await openStep(label, options, true);
    if (tx) { tx.deferAdapterCommit = true; touchTransaction(); }
    // OV-460: an edit made while another action is still running on the open
    // step applies at once and is a step under its own name: it names the
    // step from here on. A call the running action makes itself keeps the
    // name: one made synchronously inside its `fn`, or one that says so
    // (`joinsOpenStep`). Each action commits under the name it started with.
    // OV-464: a call that joins belongs to the action that started on the
    // step with its label, and commits under the name that action started
    // with, whatever a later action has named the step since.
    if (open?.actions && !actionCallDepth && !options.joinsOpenStep) open.label = label;
    const joined = options.joinsOpenStep ? open?.actions?.find((started) => started.label === label) : null;
    const action = { label, stepLabel: joined ? joined.stepLabel : tx?.label };

    try {
      const busy = mutationAvailability();
      if (busy) { cancel(tx); return busy; }
      // The step closes when an action on it ends, so the list is never cleared.
      if (tx) (tx.actions ??= []).push(action);
      const result = await callAction(fn);
      await commitIntent(tx, options, action);
      return result;
    } catch (error) {
      cancel(tx);
      throw error;
    }
  };

  /**
   * @param {string} [label]
   * @param {HistoryStepBase} [options]
   * @returns {HistoryTransaction | Promise<HistoryTransaction> | null}
   */
  const beginCheckpoint = (label = 'Artifact change', options = {}) => {
    if (restoring.value || capturing.value) return null;
    /**
     * @param {CheckpointRecord} before
     * @param {boolean | undefined} drawn
     * @returns {HistoryTransaction}
     */
    const createTransaction = (before, drawn) => {
      if (drawn !== undefined) before.drawn = drawn;
      if (currentCheckpoint !== null) setCurrentCheckpoint(before);
      emitHistoryDiagnostic({ type: 'begin', scope: 'artifact-checkpoint', label });
      return {
        label,
        before,
        // The settings a committed checkpoint changed name the live values it asks to apply.
        beforeIntent: currentIntent,
        closed: false,
        source: options.source || '',
        takeEntry: options.takeEntry,
        sessionMark: stepOpened()
      };
    };
    const captureBefore = () => {
      notifyCheckpointCapture(options, 'before-start');
      // Read as the capture starts; the snapshot service takes the figure at once.
      const drawn = options.drawn?.before();
      const before = captureCheckpoint();
      // A project step's part keeps the intent values it changes (D-15-51),
      // so it reads the intent as it begins.
      const intent = options.takeEntry ? captureIntent() : null;
      /** @param {CheckpointRecord} record @param {IntentRecord | null} intentRecord */
      const open = (record, intentRecord) => {
        if (intentRecord) setCurrentIntent(intentRecord);
        const transaction = createTransaction(record, drawn);
        notifyCheckpointCapture(options, 'before-end');
        return transaction;
      };
      return isPromiseLike(before) || isPromiseLike(intent)
        ? Promise.all([before, intent]).then(([record, intentRecord]) => open(record, intentRecord))
        : open(before, intent);
    };

    if (!openTransaction()) return captureBefore();
    return settlePendingIntent().then(captureBefore);
  };

  /**
   * @param {HistoryTransaction | null | undefined} transaction
   * @param {HistoryStepBase} [options]
   */
  const commitCheckpoint = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    notifyCheckpointCapture(options, 'after-start');
    const after = await captureCheckpoint();
    if (options.drawn) after.drawn = transaction.before.drawn === true && options.drawn.change;
    const intent = await captureIntent();
    transaction.closed = true;
    setCurrentCheckpoint(after);
    setCurrentIntent(intent);

    if (after.signature === transaction.before.signature) {
      endStep(null, transaction);
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
    const changes = buildIntentPatch(transaction.beforeIntent ?? {}, intent.intent);
    /** @type {CheckpointEntry} */
    const entry = {
      type: 'checkpoint',
      label: transaction.label || options.label || 'Artifact change',
      before: transaction.before,
      after,
      byteSize: transaction.before.byteSize + after.byteSize,
      fileIds,
      // A project step's part restores these when its drawing is not shown.
      ...(transaction.takeEntry ? { changes } : {})
    };
    endStep(entry, transaction);
    emitHistoryDiagnostic({
      type: 'commit',
      scope: 'artifact-checkpoint',
      label: entry.label,
      created: true,
      bytes: entry.byteSize
    });
    enforceLimits(entry);
    touch();
    notifyCheckpointCapture(options, 'after-end');
    settingsStepListener(changes);
    return true;
  };

  /**
   * @template R
   * @param {string} label
   * @param {() => R | Promise<R>} fn
   * @param {HistoryStepOptions<Awaited<R>>} [options]
   */
  const runUndoableCheckpoint = (label, fn, options = {}) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (typeof fn !== 'function') return undefined;
    if (restoring.value || capturing.value || activeCheckpoint) return fn();
    /** @param {HistoryTransaction | null} tx */
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
          intentGeneration += 1;
          releaseUnreferencedFiles();
          return result;
        }
        await commitCheckpoint(tx, options);
        return result;
      } catch (error) {
        if (tx) tx.closed = true;
        intentGeneration += 1;
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

  // D-24 (OV-731): a rerender drawn from the current state is shown (D-36).
  // While the last step is a checkpoint and nothing is open, the state is the
  // one its after side holds, so a side not yet drawn keeps the figure the
  // rerender drew: Redo then shows it without another Python run.
  const keepRedrawnResult = async () => {
    const entry = topEntry(undoStack);
    if (entry?.type !== 'checkpoint') return false;
    const keeps = () => entry === topEntry(undoStack) && entry.after.drawn === false
      && !openTransaction() && !activeCheckpoint && !activeReplacement && !restoring.value && !capturing.value;
    if (!keeps()) return false;
    const after = await captureCheckpoint();
    if (!keeps()) return false;
    after.drawn = true;
    removeEntryBytes(entry);
    entry.after = after;
    entry.byteSize = entry.before.byteSize + after.byteSize;
    entry.fileIds = new Set([...entry.before.fileIds, ...after.fileIds]);
    addEntryBytes(entry);
    setCurrentCheckpoint(after);
    enforceLimits(entry);
    return true;
  };

  // OV-872: what a rerender writes as its figure is shown (the Results it
  // replaces, OV-873, and the Legend rows Python drew) is drawn state, not an
  // edit: an intent step that began after the rerender was requested takes it
  // into the state it began from, so no step records it, as when no step is
  // open. Its own changes stay in it. The step that requested the rerender
  // (begun at `requestedAt` or earlier) keeps the write as its own, so its
  // Undo restores the rows from before its edit (D-24 (OV-731), D-40). Both
  // ways rest on the serial a step takes as it begins (`openStep`): a step
  // begun after the request is numbered above it, and a step whose own edit
  // made the request is numbered at or below it, also when it first waits for
  // another owner's step to settle. A request no step made, read while a step
  // is open, counts that step as the requester, so the write stays in it (as
  // before D-40; no such request source exists now). The
  // intent is read on both sides of the write in the same task, so the
  // binder's timing does not change; a snapshot service that builds the
  // intent asynchronously keeps the write in the step.
  // Arrays are leaves (`buildIntentPatch`): a drawn array replaces the whole
  // array in the step's `before`, so a write the step made to the same array
  // before the rerender binds drops out of the step (no action does so now).
  // While an intent capture runs (`capturing`, as a step begins) the write
  // is not adopted and the step that begins may record it.
  /**
   * @param {() => void} write
   * @param {{ requestedAt: number }} request `intentStepSerial()` when the rerender was requested.
   */
  const adoptDrawnWrite = (write, { requestedAt }) => {
    const open = restoring.value || capturing.value ? null : openTransaction();
    if (!open || (open.serial ?? 0) <= requestedAt) { write(); return; }
    const before = buildIntent();
    write();
    if (isPromiseLike(before)) return;
    const after = buildIntent();
    diagnostics.intentBuilds += 2;
    if (isPromiseLike(after)) return;
    const changes = buildIntentPatch(before, after);
    if (changes.length === 0) return;
    open.before = applyIntentPatch(open.before, changes, 'redo');
    open.beforeSignature = computeSignature(open.before, 'intent');
  };

  /**
   * @param {HistoryArtifactHandle} handle
   * @param {'before' | 'after'} phase
   */
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

  /** @param {'before' | 'after'} phase */
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

  /**
   * @param {HistoryStepBase | undefined} options
   * @param {'before' | 'after'} phase
   * @returns {Promise<ReplacementIntentCheckpoint>}
   */
  const captureReplacementIntentCheckpoint = async (options, phase) => {
    if (typeof options?.captureIntentCheckpoint !== 'function') {
      return { enabled: false, value: undefined };
    }
    return {
      enabled: true,
      value: await options.captureIntentCheckpoint({ phase })
    };
  };

  /**
   * @param {ReplacementIntentCheckpoint | undefined} checkpoint
   * @param {HistoryStepBase['restoreIntentCheckpoint'] | null | undefined} restore
   */
  const restoreReplacementIntentCheckpoint = async (checkpoint, restore) => {
    if (!checkpoint?.enabled || typeof restore !== 'function') return;
    restoring.value = true;
    try {
      await restore(checkpoint.value);
    } finally {
      restoring.value = false;
    }
  };

  /**
   * @param {string} [label]
   * @param {HistoryStepBase} [options]
   * @returns {Promise<HistoryTransaction | null>}
   */
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

  /**
   * @param {HistoryArtifactHandle} before
   * @param {HistoryArtifactHandle} after
   */
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

  /**
   * @param {HistoryArtifactHandle | undefined} handle
   * @returns {Set<string>}
   */
  const artifactHandleFileIds = (handle) => new Set(handle?.fileIds || []);

  /**
   * @param {HistoryTransaction | null | undefined} transaction
   * @param {HistoryStepBase} [options]
   */
  const commitAppliedArtifactReplacement = async (transaction, options = {}) => {
    if (!transaction || transaction.closed) return false;
    notifyCheckpointCapture(options, 'after-start');
    const after = await captureArtifactHandle('after');
    const afterIntentCheckpoint = await captureReplacementIntentCheckpoint(
      options,
      'after'
    );
    transaction.closed = true;
    intentGeneration += 1;
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
    /** @type {ArtifactReplacementEntry} */
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
    endStep(entry, transaction);
    enforceLimits(entry);
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

  /**
   * @param {HistoryArtifactHandle} handle
   * @param {{ refreshIntent?: boolean }} [options]
   */
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

  /**
   * @template R
   * @param {string} label
   * @param {(before: HistoryArtifactHandle | null) => R | Promise<R>} fn
   * @param {HistoryStepOptions<Awaited<R>>} [options]
   */
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
        intentGeneration += 1;
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

  /**
   * @param {string | undefined} label
   * @param {UndoableCommand | null | undefined} command
   * @returns {CommandEntry | null}
   */
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

  /** @param {unknown} result */
  const commandSucceeded = (result) => result !== false;

  /**
   * @param {{ label?: string } | null | undefined} entry
   * @param {string} action
   */
  const warnCommandNotApplied = (entry, action) => {
    const label = entry?.label ? ` "${entry.label}"` : '';
    console.warn(`Undo/redo command${label} could not ${action}; history was left unchanged.`);
  };

  const refreshCurrentIntent = async () => {
    setCurrentIntent(await captureIntent());
  };

  /**
   * @param {string} label
   * @param {() => UndoableCommand | null | undefined | Promise<UndoableCommand | null | undefined>} buildCommand
   */
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
    endStep(command, command);
    enforceLimits(command);
    touch();
    return true;
  };

  /**
   * @param {IntentEntry} entry
   * @param {'undo' | 'redo'} direction
   */
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
    // The restored live setting values return to the Result once applied (6.4).
    settingsStepListener(entry.changes);
  };

  /** @param {CheckpointRecord} record */
  const applyCheckpointEntry = async (record) => {
    restoring.value = true;
    try {
      await applyCheckpoint(record.checkpoint, { drawn: record.drawn === true });
    } finally {
      restoring.value = false;
    }
    // A checkpoint's figure can predate the apply of its live values (the
    // Redo of Reset Settings); one already equal applies nothing.
    settingsStepListener(null);
    setCurrentCheckpoint(record);
    await refreshCurrentIntent();
  };

  /**
   * @param {ArtifactReplacementEntry} entry
   * @param {'undo' | 'redo'} direction
   */
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

  /**
   * @param {CommandEntry} entry
   * @param {'undo' | 'redo'} direction
   */
  const applyCommandWithFlag = async (entry, direction) => {
    restoring.value = true;
    try {
      const result = direction === 'undo' ? await entry.revert() : await entry.apply();
      return commandSucceeded(result);
    } finally {
      restoring.value = false;
    }
  };

  // D-15-53: the Session-level values the step changed with this drawing's.
  /** @param {HistoryEntry} entry @param {'undo' | 'redo'} direction */
  const restoreSessionPart = async (entry, direction) => {
    if (!entry.session) return;
    restoring.value = true;
    try {
      await applySessionPart(entry.session, direction);
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

  // OV-1339: the busy answer for the step an Undo (or Redo) takes next, the
  // top of its stack. The action and its enabled state read this one check.
  /** @param {HistoryEntry[]} stack */
  const topStepAvailability = (stack) => (
    stack.length > 0 ? stepAvailability(stack[stack.length - 1]) : null
  );

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
    const stepBusy = topStepAvailability(undoStack);
    if (stepBusy) return stepBusy;
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
      await applyArtifactReplacementEntry(entry, 'undo');
    } else {
      await applyIntentEntry(entry, 'undo');
    }
    await restoreSessionPart(entry, 'undo');
    undoStack.pop();
    redoStack.push(entry);
    enforceLimits(entry);
    touch();
    return true;
  });

  const redo = () => traverse(async () => {
    const busy = await settleBeforeTraversal();
    if (busy) return busy;
    if (restoring.value || redoStack.length === 0) return false;
    const stepBusy = topStepAvailability(redoStack);
    if (stepBusy) return stepBusy;
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
      await applyArtifactReplacementEntry(entry, 'redo');
    } else {
      await applyIntentEntry(entry, 'redo');
    }
    await restoreSessionPart(entry, 'redo');
    redoStack.pop();
    undoStack.push(entry);
    enforceLimits(entry);
    touch();
    return true;
  });

  // D-15-51: this drawing's part of a project step (Reset all drawings). An
  // intent part restores only the values still as the step left them (a
  // later step of this drawing on the same path is kept, OIC-030); a value
  // a redraw only laid out anew is as left (`editedPartOf`, OV-1655). A
  // checkpoint part restores its checkpoint when `whole` (its drawing is
  // shown and kept every newer step, OV-1169), and otherwise the intent
  // values it changed; the Result shows them when the drawing is shown again
  // (departed signature).
  /**
   * @param {IntentEntry | CheckpointEntry} entry
   * @param {'undo' | 'redo'} direction
   * @param {boolean} whole
   */
  const applyProjectPart = (entry, direction, whole) => traverse(async () => {
    if (entry.type === 'checkpoint' && whole) {
      await applyCheckpointEntry(direction === 'undo' ? entry.before : entry.after);
    } else {
      await refreshCurrentIntent();
      const changes = changesAsLeft(currentIntent, entry.changes || [], direction, editedPartOf);
      if (changes.length > 0) await applyIntentEntry({ ...entry, type: 'intent', changes }, direction);
    }
    touch();
  });
  // The router's own traversal steps (a project step, a Session-level value)
  // run while this History restores, so no edit they cause records a step.
  /** @param {() => unknown} fn */
  const restoreWith = (fn) => traverse(async () => {
    restoring.value = true;
    try { await fn(); } finally { restoring.value = false; }
    touch();
  });
  /** @param {'undo' | 'redo'} stack */
  const topSequence = (stack) => topEntry(stack === 'undo' ? undoStack : redoStack)?.sequence ?? null;
  const dropRedo = () => {
    if (redoStack.length === 0) return;
    clearRedo();
    touch();
  };

  /** @param {HistoryEntry[]} stack */
  const canTraverse = (stack) => (
    stack.length > 0 && !restoring.value && !historyAvailability() && !topStepAvailability(stack)
  );
  const canUndo = () => canTraverse(undoStack);
  const canRedo = () => canTraverse(redoStack);
  const undoLabel = () => (canUndo() ? undoStack[undoStack.length - 1].label : '');
  const redoLabel = () => (canRedo() ? redoStack[redoStack.length - 1].label : '');

  /**
   * The one listener of the settings a History step records, a checkpoint
   * commits, or an Undo or Redo restores (`null`: any setting, a checkpoint
   * restore). The composition root registers it once both History and the
   * automatic rerender exist (R1 (c), PD-OI-102).
   * @param {(changes: readonly HistoryIntentChange[] | null) => void} listener
   */
  const setSettingsStepListener = (listener) => { settingsStepListener = listener; };

  return {
    setSettingsStepListener,
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
    markChangeOutsideHistory,
    historyAvailability,
    keepRedrawnResult,
    adoptDrawnWrite,
    intentStepSerial: () => intentStepSerial,
    traversalPending: () => traversals.value > 0,
    initializeIntentBaseline,
    markSavePoint,
    redo,
    redoLabel,
    applyProjectPart,
    restoreWith,
    topSequence,
    lostStepAfter: (/** @type {number} */ sequence) => droppedThrough > sequence,
    dropRedo,
    settleBeforeTraversal,
    // PD-OI-090: a drawing that is shown again takes the shared selections in
    // its intent that changed meanwhile as its current intent, so its Undo and
    // Redo restore only what its own steps changed.
    refreshCurrentIntent,
    restoring,
    capturing,
    revision,
    runUndoable,
    runUndoableArtifactReplacement,
    runUndoableCheckpoint,
    runUndoableCommand,
    // R11: the open step ends now, in this History (a drawing switch).
    settlePendingIntent,
    // The Session values the open step opened with, if a step is open.
    openStepMark: () => openTransaction()?.sessionMark ?? null,
    heldIdentities,
    undo,
    undoLabel
  };
};

/**
 * @typedef {{ drawingId: string, entry: IntentEntry | CheckpointEntry }} ProjectPart
 * @typedef {{
 *   type: 'project', label: string, sequence: number, parts: ProjectPart[], session: IntentEntry | null,
 *   byteSize: number, fileIds: Set<string>
 * }} ProjectEntry
 *   A step of the project stream: its part in each drawing it changed (Reset
 *   all drawings) and the Session-level values it changed (R15-23).
 * @typedef {object} SessionHistoryPorts
 *   The capture and restore of the Session-level values, which no drawing owns.
 * @property {() => HistoryIntent} buildIntent Synchronous: a step's end reads it.
 * @property {(intent: HistoryIntent, context: { changes: HistoryIntentChange[], direction: 'undo' | 'redo' }) => unknown} applyIntent
 * @typedef {HistoryBudgetOptions & {
 *   drawingIds: readonly string[],
 *   activeDrawingId: () => string,
 *   managerOptions: (drawingId: string) => HistoryManagerOptions,
 *   sessionHistory: SessionHistoryPorts,
 *   makeRef?: <T>(value: T) => { value: T },
 *   maxActions?: number,
 *   onHistoryLimit?: (message: string) => void
 * }} HistoryRouterOptions
 *   The drawings, each with its own History; the shown drawing, whose History
 *   Undo and Redo walk; the capture and restore of one drawing's History and
 *   of the Session-level values; and the budget every History shares.
 */

// The budget member of the project stream (D-15-52).
const PROJECT_STREAM = Symbol('project stream');

// PD-OI-090: each drawing has its own History. The router keeps the manager
// API and routes it to the shown drawing's manager; a step that a call opened
// closes in the manager that opened it, also after a drawing switch. A save
// point and the baselines cover every drawing; the settings-step listener
// hears the shown drawing's steps; a mutation or traversal pending in any
// drawing is pending. The managers share one byte budget and one file
// retention (`createHistoryBudget`).
//
// D-15-51, D-15-52, R15-23, D-15-53: the router also keeps the project
// stream. A step of it is a project operation (Reset all drawings: a part in
// each drawing) or a step that changed Session-level values alone (the
// Session title, the app settings, the LOSAT execution, the drafts). The
// router reads those values at the end of every step: a step that changed a
// drawing value is that drawing's step whole, and the Session-level values it
// changed are its entry's `session` part, restored with it where still as it
// left them. Undo takes the newer of the shown drawing's and the project's
// top steps, Redo the older. A new project step clears every Redo; a new
// drawing step clears its own and the project's. Eviction drops hidden
// drawings' steps first, then project steps, then the shown drawing's.
/** @param {HistoryRouterOptions} options */
export const createHistoryRouter = ({
  drawingIds,
  activeDrawingId,
  managerOptions,
  sessionHistory,
  makeRef = makeBox,
  maxActions = DEFAULT_MAX_ACTIONS,
  onHistoryLimit = () => {},
  ...budgetOptions
}) => {
  /** @type {Map<HistoryBudgetMember, string | typeof PROJECT_STREAM>} The stream of each budget member. */
  const memberStreams = new Map();
  /** @param {HistoryBudgetMember} member */
  const evictionRank = (member) => {
    const stream = memberStreams.get(member);
    if (stream === PROJECT_STREAM) return 1;
    return stream === activeDrawingId() ? 2 : 0;
  };
  const budget = createHistoryBudget({
    ...budgetOptions,
    evictionOrder: (members) => [...members].sort((left, right) => evictionRank(left) - evictionRank(right))
  });
  /** @param {string | typeof PROJECT_STREAM} stream @param {HistoryBudgetMember} member */
  const joinBudget = (stream, member) => { memberStreams.set(member, stream); budget.join(member); };

  let sequence = 0;
  /** @type {ProjectEntry[]} */
  const projectUndo = [];
  /** @type {ProjectEntry[]} */
  const projectRedo = [];
  const projectRevision = makeRef(0);
  const touchProject = () => { projectRevision.value += 1; };
  /** @returns {SessionMark} */
  const sessionRecord = () => {
    const intent = cloneIntentValue(sessionHistory.buildIntent());
    return { intent, signature: JSON.stringify(intent) };
  };
  /** @param {SessionMark} record @param {HistoryIntentChange[]} changes @param {'undo' | 'redo'} [direction] */
  const patchedRecord = (record, changes, direction = 'redo') => {
    const intent = applyIntentPatch(record.intent, changes, direction);
    return { intent, signature: JSON.stringify(intent) };
  };
  /** @param {SessionMark} from @param {SessionMark} to */
  const recordChanges = (from, to) => (from.signature === to.signature ? [] : buildIntentPatch(from.intent, to.intent));
  // The Session-level values as the last step left them.
  let sessionCurrent = sessionRecord();

  /** @param {ProjectEntry} entry */
  const pushProject = (entry) => {
    projectUndo.push(entry);
    projectRedo.length = 0;
    every().forEach((manager) => manager.dropRedo());
    // The step-count window drops silently (OV-1365).
    while (projectUndo.length > maxActions) projectUndo.shift();
    budget.enforce(entry);
    touchProject();
  };
  // The Session-level values a step changed. A step that changed a drawing
  // value (`mark` null) takes every change since the last step ended, so it
  // keeps one made outside any step (D-15-53 2). A step that changed none
  // takes only the changes since it opened (`mark`); an earlier one stays
  // for the next step.
  /** @param {string} label @param {SessionMark | null} [mark] @returns {IntentEntry | null} */
  const sessionStep = (label, mark = null) => {
    const after = sessionRecord();
    const changes = recordChanges(mark ?? sessionCurrent, after);
    sessionCurrent = mark ? patchedRecord(sessionCurrent, changes) : after;
    return changes.length > 0
      ? { type: 'intent', label, changes, byteSize: JSON.stringify(changes).length * 2, fileIds: new Set() }
      : null;
  };
  // The project step of a step's end: its parts and the Session-level values
  // it changed; none when it has neither.
  /** @param {string} label @param {number} stamp @param {ProjectPart[]} parts @param {IntentEntry | null} session */
  const recordProjectStep = (label, stamp, parts, session) => {
    if (!session && parts.length === 0) return;
    /** @type {Set<string>} */
    const fileIds = new Set();
    parts.forEach(({ entry }) => entry.fileIds.forEach((id) => fileIds.add(id)));
    pushProject({
      type: 'project', label, sequence: stamp, parts, session, fileIds,
      byteSize: parts.reduce((bytes, { entry }) => bytes + entry.byteSize, session?.byteSize ?? 0)
    });
  };
  /** @param {HistoryEntry | null} entry @param {string} label @param {SessionMark | null} mark @returns {IntentEntry | null} */
  const stepEnded = (entry, label, mark) => {
    const stamp = ++sequence;
    const session = sessionStep(label, entry ? null : mark);
    if (!entry) {
      recordProjectStep(label, stamp, [], session);
      return null;
    }
    entry.sequence = stamp;
    if (projectRedo.length > 0) { projectRedo.length = 0; touchProject(); }
    return session;
  };

  const managers = new Map(drawingIds.map((drawingId) => /** @type {const} */ ([drawingId, createHistoryManager({
    ...managerOptions(drawingId),
    onHistoryLimit,
    stepOpened: sessionRecord,
    stepEnded,
    applySessionPart: (part, direction) => applySession(part, direction, { asLeft: true }),
    budget: {
      ...budget,
      join: (member) => joinBudget(drawingId, member),
      // OV-1260 keeps the step the user just recorded, undid or redid in the
      // shown drawing; a hidden drawing's step drops in the order that follows
      // the shown drawing (A3 review Low 4).
      enforce: (kept) => budget.enforce(drawingId === activeDrawingId() ? kept : null)
    }
  })])));
  /** @param {string} drawingId */
  const managerOf = (drawingId) => {
    const manager = managers.get(drawingId);
    if (!manager) throw new Error(`No History exists for drawing ${drawingId}.`);
    return manager;
  };
  const every = () => [...managers.values()];
  const others = () => [...managers].filter(([drawingId]) => drawingId !== activeDrawingId()).map(([, manager]) => manager);
  const template = managerOf(drawingIds[0]);
  /** @typedef {typeof template} DrawingHistory */
  const projectEntries = () => [...projectUndo, ...projectRedo];
  joinBudget(PROJECT_STREAM, {
    entryBytes: () => projectEntries().reduce((bytes, entry) => bytes + entry.byteSize, 0),
    retainedFileIds: () => {
      /** @type {Set<string>} */
      const ids = new Set();
      projectEntries().forEach((entry) => entry.fileIds.forEach((id) => ids.add(id)));
      return ids;
    },
    dropOldestStep: (kept) => {
      if (projectUndo.length > 0 && projectUndo[0] !== kept) { projectUndo.shift(); touchProject(); return 'undo'; }
      const keptRedo = projectRedo.filter((entry) => entry === kept);
      if (keptRedo.length === projectRedo.length) return null;
      projectRedo.splice(0, projectRedo.length, ...keptRedo);
      touchProject();
      return 'redo';
    },
    noteEvicted: onHistoryLimit
  });

  // A project step's Session-level values, or a drawing entry's (`asLeft`:
  // only where still as its step left them). A value changed outside any
  // step stays a change for the next step (D-15-53).
  /** @param {IntentEntry} entry @param {'undo' | 'redo'} direction @param {{ asLeft?: boolean }} [options] */
  const applySession = async (entry, direction, { asLeft = false } = {}) => {
    const live = sessionRecord().intent;
    const changes = asLeft ? changesAsLeft(live, entry.changes, direction) : entry.changes;
    if (changes.length === 0) return;
    await sessionHistory.applyIntent(applyIntentPatch(live, changes, direction), { changes, direction });
    sessionCurrent = patchedRecord(sessionCurrent, changes, direction);
  };
  // Every drawing settles its open step and none restores; the Session-level
  // values apply first (a drawing's queued palette restores under the
  // Instant Preview flag it was queued with, OV-1160), then the parts of the
  // drawings that are not shown, the shown drawing's last, while the shown
  // drawing's History restores (no edit they cause records a step).
  /** @param {ProjectEntry} entry @param {'undo' | 'redo'} direction */
  const traverseProject = async (entry, direction) => {
    for (const manager of every()) {
      const busy = await manager.settleBeforeTraversal();
      if (busy) return busy;
    }
    const [from, to] = direction === 'undo' ? [projectUndo, projectRedo] : [projectRedo, projectUndo];
    if (every().some((manager) => manager.restoring.value) || from[from.length - 1] !== entry) return false;
    const shownId = activeDrawingId();
    const shown = managerOf(shownId);
    await shown.restoreWith(async () => {
      if (entry.session) await applySession(entry.session, direction);
      for (const part of entry.parts) {
        if (part.drawingId !== shownId) await managerOf(part.drawingId).applyProjectPart(part.entry, direction, false);
      }
      const own = entry.parts.find((part) => part.drawingId === shownId);
      if (own) await shown.applyProjectPart(own.entry, direction, !shown.lostStepAfter(entry.sequence));
    });
    from.pop();
    to.push(entry);
    budget.enforce(entry);
    touchProject();
    return true;
  };
  // OV-1339: an Undo or Redo is pending from its call, before the shown
  // drawing settles its open step, so the drawing switch waits for it.
  const projectTraversals = makeRef(0);
  /** @param {'undo' | 'redo'} direction */
  const traverse = async (direction) => {
    projectTraversals.value += 1;
    try {
      return await traverseShown(direction);
    } finally {
      projectTraversals.value -= 1;
    }
  };
  /** @param {'undo' | 'redo'} direction */
  const traverseShown = async (direction) => {
    const shown = managerOf(activeDrawingId());
    const busy = await shown.settleBeforeTraversal();
    if (busy) return busy;
    const stack = direction === 'undo' ? projectUndo : projectRedo;
    const project = stack[stack.length - 1];
    const own = shown.topSequence(direction);
    const ownFirst = (/** @type {number} */ ownSequence) => (
      direction === 'undo' ? ownSequence > project.sequence : ownSequence < project.sequence);
    if (!project || (own !== null && ownFirst(own))) return direction === 'undo' ? shown.undo() : shown.redo();
    return traverseProject(project, direction);
  };
  /** @param {'undo' | 'redo'} direction The project step that Undo or Redo takes next, if any. */
  const nextProjectStep = (direction) => {
    const stack = direction === 'undo' ? projectUndo : projectRedo;
    const project = stack[stack.length - 1];
    if (!project) return null;
    const own = managerOf(activeDrawingId()).topSequence(direction);
    if (own === null) return project;
    return (direction === 'undo' ? own > project.sequence : own < project.sequence) ? null : project;
  };
  const shownIdle = () => {
    const shown = managerOf(activeDrawingId());
    return !shown.restoring.value && !shown.historyAvailability();
  };
  /** @param {string} drawingId */
  const projectPartsOf = (drawingId) => projectEntries()
    .flatMap((entry) => entry.parts.filter((part) => part.drawingId === drawingId).map((part) => part.entry));
  /** @type {{ undo: ProjectEntry | null, redo: ProjectEntry | null }} */
  let projectSavePoint = { undo: null, redo: null };
  const resetProjectStream = () => {
    projectUndo.length = 0;
    projectRedo.length = 0;
    sessionCurrent = sessionRecord();
    projectSavePoint = { undo: null, redo: null };
    touchProject();
  };

  const own = {
    managerOf,
    revision: { get value() { return projectRevision.value + managerOf(activeDrawingId()).revision.value; } },
    undo: () => traverse('undo'),
    redo: () => traverse('redo'),
    canUndo: () => managerOf(activeDrawingId()).canUndo() || (projectUndo.length > 0 && shownIdle()),
    canRedo: () => managerOf(activeDrawingId()).canRedo() || (projectRedo.length > 0 && shownIdle()),
    undoLabel: () => (shownIdle() ? nextProjectStep('undo')?.label ?? managerOf(activeDrawingId()).undoLabel() : ''),
    redoLabel: () => (shownIdle() ? nextProjectStep('redo')?.label ?? managerOf(activeDrawingId()).redoLabel() : ''),
    getUndoCount: () => managerOf(activeDrawingId()).getUndoCount() + projectUndo.length,
    getRedoCount: () => managerOf(activeDrawingId()).getRedoCount() + projectRedo.length,
    // A drawing that is not shown has no artifact installed, so its baseline
    // is its intent only.
    captureBaseline: async (/** @type {string} */ label = 'Baseline') => {
      await managerOf(activeDrawingId()).captureBaseline(label);
      for (const manager of others()) await manager.initializeIntentBaseline(label);
      resetProjectStream();
    },
    initializeIntentBaseline: async (
      /** @type {string} */ label = 'Intent baseline',
      /** @type {{ isCurrent?: () => boolean }} */ options = {}
    ) => {
      for (const manager of every()) {
        if (!await manager.initializeIntentBaseline(label, options)) return false;
      }
      resetProjectStream();
      return true;
    },
    // D-15-51: an action that changes every drawing (Reset all drawings) is
    // one project step. `fn` runs in this task as the shown drawing's
    // checkpoint, as a caller of the action expects; each other drawing
    // captured its state first and records what `fn` changed in it. Their
    // entries are the step's parts, with the Session-level values it changed.
    // `shouldCommit` answers, from the result, whether the step is recorded.
    /**
     * @template R
     * @param {string} label
     * @param {() => R | Promise<R>} fn
     * @param {{ shouldCommit?: (result: unknown) => boolean }} [options]
     */
    runProjectStep: (label, fn, { shouldCommit = () => true } = {}) => {
      const shownId = activeDrawingId();
      /** @type {ProjectPart[]} */
      const parts = [];
      /** @param {string} drawingId @returns {(entry: HistoryEntry | null) => void} */
      const take = (drawingId) => (entry) => {
        if (entry?.type === 'intent' || entry?.type === 'checkpoint') parts.push({ drawingId, entry });
      };
      const hidden = [...managers].filter(([drawingId]) => drawingId !== shownId);
      const opening = hidden.map(([drawingId, manager]) => manager.begin(label, { source: 'project', takeEntry: take(drawingId) }));
      /** @param {boolean} applied */
      const settle = async (applied) => {
        const steps = await Promise.all(opening);
        for (const [index, [, manager]] of hidden.entries()) {
          if (applied) await manager.commit(steps[index]); else manager.cancel(steps[index]);
        }
        if (applied) recordProjectStep(label, ++sequence, parts, sessionStep(label));
      };
      /** @type {ReturnType<DrawingHistory['runUndoableCheckpoint']>} */
      let result;
      try {
        result = managerOf(shownId).runUndoableCheckpoint(label, fn, { takeEntry: take(shownId), shouldCommit });
      } catch (error) {
        void settle(false);
        throw error;
      }
      if (!isPromiseLike(result)) {
        void settle(shouldCommit(result));
        return result;
      }
      return Promise.resolve(result).then(
        async (value) => { await settle(shouldCommit(value)); return value; },
        async (error) => { await settle(false); throw error; }
      );
    },
    markSavePoint: () => {
      every().forEach((manager) => manager.markSavePoint());
      projectSavePoint = { undo: projectUndo[projectUndo.length - 1] ?? null, redo: projectRedo[projectRedo.length - 1] ?? null };
    },
    // D-15-53: a drawing switch resets Session-level transient values (the
    // feature panel tab) without a step. Only what `reset` changes becomes
    // recorded, as the last step left it and as an open step opened with it;
    // an open step's own edits stay its changes.
    adoptSessionReset: (/** @type {() => void} */ reset) => {
      const before = sessionRecord();
      reset();
      const changes = recordChanges(before, sessionRecord());
      if (changes.length === 0) return;
      sessionCurrent = patchedRecord(sessionCurrent, changes);
      every().forEach((manager) => {
        const mark = manager.openStepMark();
        if (mark) Object.assign(mark, patchedRecord(mark, changes));
      });
    },
    // OV-110: the Result identities any drawing's History holds.
    heldIdentities: () => {
      /** @type {Set<string>} */
      const ids = new Set();
      managers.forEach((manager, drawingId) => manager.heldIdentities(projectPartsOf(drawingId)).forEach((id) => ids.add(id)));
      return ids;
    },
    // D-15-53: a Session-level value changed outside any step is no change
    // until a step takes it; one an open step changed since it opened is (as
    // a manager counts its own values).
    hasChangesSinceSavePoint: async () => {
      if ((projectUndo[projectUndo.length - 1] ?? null) !== projectSavePoint.undo) return true;
      if ((projectRedo[projectRedo.length - 1] ?? null) !== projectSavePoint.redo) return true;
      const marks = every().flatMap((manager) => manager.openStepMark() ?? []);
      if (marks.length > 0) {
        const { signature } = sessionRecord();
        if (marks.some((mark) => mark.signature !== signature)) return true;
      }
      for (const manager of every()) if (await manager.hasChangesSinceSavePoint()) return true;
      return false;
    },
    /** @param {Parameters<DrawingHistory['setSettingsStepListener']>[0]} listener */
    setSettingsStepListener: (listener) => managers.forEach((manager, drawingId) => manager.setSettingsStepListener(
      (changes) => { if (drawingId === activeDrawingId()) listener(changes); }
    )),
    mutationPending: () => every().some((manager) => manager.mutationPending()),
    traversalPending: () => projectTraversals.value > 0 || every().some((manager) => manager.traversalPending())
  };

  /** @type {WeakMap<object, DrawingHistory>} The manager that opened each step. */
  const openedBy = new WeakMap();
  /** @param {DrawingHistory} manager @param {string} key */
  const memberOf = (manager, key) => /** @type {(...args: unknown[]) => unknown} */ (
    /** @type {Record<string, unknown>} */ (/** @type {unknown} */ (manager))[key]);
  /** @param {string} key */
  const activeMember = (key) => (
    /** @type {Record<string, unknown>} */ (/** @type {unknown} */ (managerOf(activeDrawingId()))))[key];
  const OPENERS = new Set(['begin', 'beginCheckpoint', 'beginArtifactReplacement']);
  const CLOSERS = new Set(['commit', 'cancel', 'commitCheckpoint', 'commitAppliedArtifactReplacement']);
  /** @type {Record<string, unknown>} */
  const routed = {};
  Object.entries(template).forEach(([key, value]) => {
    if (Object.hasOwn(own, key)) return;
    if (OPENERS.has(key)) {
      routed[key] = (/** @type {unknown[]} */ ...args) => {
        const manager = managerOf(activeDrawingId());
        const opened = memberOf(manager, key)(...args);
        /** @param {unknown} step */
        const note = (step) => { if (step && typeof step === 'object') openedBy.set(step, manager); };
        if (isPromiseLike(opened)) opened.then(note, () => {}); else note(opened);
        return opened;
      };
    } else if (CLOSERS.has(key)) {
      routed[key] = (/** @type {unknown} */ step, /** @type {unknown[]} */ ...rest) => {
        const opener = step && typeof step === 'object' ? openedBy.get(step) : undefined;
        return memberOf(opener ?? managerOf(activeDrawingId()), key)(step, ...rest);
      };
    } else if (typeof value === 'function') {
      routed[key] = (/** @type {unknown[]} */ ...args) => (
        /** @type {(...args: unknown[]) => unknown} */ (activeMember(key)))(...args);
    } else {
      Object.defineProperty(routed, key, { enumerable: true, get: () => activeMember(key) });
    }
  });
  return Object.assign(/** @type {Omit<DrawingHistory, keyof typeof own>} */ (/** @type {unknown} */ (routed)), own);
};
