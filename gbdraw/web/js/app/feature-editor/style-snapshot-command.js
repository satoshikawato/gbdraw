const isPlainObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const isPromiseLike = (value) => Boolean(value && typeof value.then === 'function');

const stableValue = (value) => {
  if (Array.isArray(value)) return value.map(stableValue);
  if (!isPlainObject(value)) return value;
  return Object.fromEntries(
    Object.keys(value).sort().map((key) => [key, stableValue(value[key])])
  );
};

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const utf8Length = (value) => {
  const source = String(value ?? '');
  if (typeof TextEncoder === 'function') return new TextEncoder().encode(source).byteLength;
  // TextEncoder is available in supported browsers. This exact fallback handles
  // scalar values without allocating a second encoded copy.
  let bytes = 0;
  for (const character of source) {
    const code = character.codePointAt(0);
    bytes += code <= 0x7f ? 1 : (code <= 0x7ff ? 2 : (code <= 0xffff ? 3 : 4));
  }
  return bytes;
};

export const exactJsonByteLength = (value) => (
  utf8Length(JSON.stringify(stableValue(value === undefined ? null : value)))
);

export const immutableStyleCommandSnapshot = (value) => deepFreeze(cloneJson(value));

const preparedRetainedBytes = (prepared) => {
  if (!prepared) return 0;
  if (Number.isFinite(Number(prepared.retainedBytes)) && Number(prepared.retainedBytes) >= 0) {
    return Number(prepared.retainedBytes);
  }
  if (Object.prototype.hasOwnProperty.call(prepared, 'retainedForHistory')) {
    return exactJsonByteLength(prepared.retainedForHistory);
  }
  return 0;
};

const commandNoop = ({ label, metadata, counters }) => Object.freeze({
  label,
  metadata: deepFreeze(cloneJson(metadata || {})),
  counters: deepFreeze(cloneJson(counters || {})),
  noop: true,
  apply: () => false,
  revert: () => false,
  compensate: () => true,
  snapshot: () => null,
  estimateBytes: () => 0
});

/**
 * Build one History-compatible command around immutable semantic sides.
 *
 * `prepareTransition` may be asynchronous and must be mutation-free. The
 * command validates again after preparation. `commitTransition`,
 * `captureExactState`, and `restoreExactState` must be synchronous so there is
 * no observable await inside the commit/rollback boundary.
 */
export const buildStyleSnapshotCommand = async ({
  label = 'Change feature style',
  before,
  after,
  metadata = {},
  counters = {},
  isNoop = null,
  validateCurrent,
  prepareTransition,
  commitTransition,
  captureExactState,
  restoreExactState,
  assertCommitted = null,
  estimateRetainedBytes = null
} = {}) => {
  if (typeof validateCurrent !== 'function') {
    throw new Error('A style command requires current-state validation.');
  }
  if (typeof prepareTransition !== 'function') {
    throw new Error('A style command requires mutation-free transition preparation.');
  }
  if (typeof commitTransition !== 'function') {
    throw new Error('A style command requires a synchronous commit owner.');
  }
  if (typeof captureExactState !== 'function' || typeof restoreExactState !== 'function') {
    throw new Error('A style command requires exact capture and restore owners.');
  }

  const immutableBefore = immutableStyleCommandSnapshot(before);
  const immutableAfter = immutableStyleCommandSnapshot(after);
  const immutableMetadata = immutableStyleCommandSnapshot(metadata || {});
  const immutableBaseCounters = immutableStyleCommandSnapshot(counters || {});
  const noop = typeof isNoop === 'function'
    ? Boolean(isNoop(immutableBefore, immutableAfter))
    : JSON.stringify(stableValue(immutableBefore)) === JSON.stringify(stableValue(immutableAfter));
  if (noop) return commandNoop({ label, metadata: immutableMetadata, counters: immutableBaseCounters });

  const sideFor = (direction) => (
    direction === 'undo'
      ? { from: immutableAfter, to: immutableBefore }
      : { from: immutableBefore, to: immutableAfter }
  );

  const prepare = async (direction) => {
    const sides = sideFor(direction);
    const firstValidation = validateCurrent({ ...sides, direction, phase: 'before-prepare' });
    if (isPromiseLike(firstValidation)) {
      throw new Error('Style command validation must be synchronous.');
    }
    if (firstValidation === false) return null;
    const prepared = await prepareTransition({ ...sides, direction });
    const secondValidation = validateCurrent({
      ...sides,
      direction,
      phase: 'after-prepare',
      prepared
    });
    if (isPromiseLike(secondValidation)) {
      throw new Error('Style command validation must be synchronous.');
    }
    if (secondValidation === false) return null;
    return prepared;
  };

  // The first forward transition is fully prepared before History can inspect
  // or apply the command. This also makes its retained byte estimate fixed.
  let preparedForward = await prepare('apply');
  if (!preparedForward) {
    throw new Error('The Feature style command became stale during preflight.');
  }
  const immutableCounters = immutableStyleCommandSnapshot({
    ...immutableBaseCounters,
    ...(isPlainObject(preparedForward.counters) ? preparedForward.counters : {})
  });
  const retainedPreparedArtifacts = preparedForward?.retainedForHistory || null;

  const retainedEnvelope = {
    before: immutableBefore,
    after: immutableAfter,
    metadata: immutableMetadata
  };
  const semanticBytes = exactJsonByteLength(retainedEnvelope);
  const artifactBytes = typeof estimateRetainedBytes === 'function'
    ? Number(estimateRetainedBytes({
        before: immutableBefore,
        after: immutableAfter,
        prepared: preparedForward,
        metadata: immutableMetadata
      }))
    : preparedRetainedBytes(preparedForward);
  if (!Number.isFinite(artifactBytes) || artifactBytes < 0) {
    throw new Error('Style command retained byte estimate must be a non-negative finite number.');
  }
  const estimatedBytes = semanticBytes + artifactBytes;

  const restore = (snapshot, context) => {
    const result = restoreExactState(snapshot, context);
    if (isPromiseLike(result)) {
      throw new Error('Style command rollback must be synchronous.');
    }
    if (result === false) {
      throw new Error('Style command exact rollback was rejected.');
    }
    return true;
  };

  const execute = async (direction) => {
    // History keeps apply/revert, so this reference keeps the accounted exact
    // Result/admission envelope alive after command normalization.
    void retainedPreparedArtifacts;
    const sides = sideFor(direction);
    let prepared;
    if (direction !== 'undo' && preparedForward) {
      const valid = validateCurrent({
        ...sides,
        direction,
        phase: 'before-commit',
        prepared: preparedForward
      });
      if (isPromiseLike(valid)) throw new Error('Style command validation must be synchronous.');
      if (valid === false) return false;
      prepared = preparedForward;
      preparedForward = null;
    } else {
      prepared = await prepare(direction);
      if (!prepared) return false;
      const valid = validateCurrent({
        ...sides,
        direction,
        phase: 'before-commit',
        prepared
      });
      if (isPromiseLike(valid)) throw new Error('Style command validation must be synchronous.');
      if (valid === false) return false;
    }

    const snapshot = captureExactState({ direction, prepared });
    if (isPromiseLike(snapshot)) {
      throw new Error('Style command capture must be synchronous.');
    }
    try {
      const committed = commitTransition({ ...sides, direction, prepared });
      if (isPromiseLike(committed)) {
        throw new Error('Style command commit must be synchronous.');
      }
      if (committed === false) {
        restore(snapshot, { direction, prepared, reason: 'rejected-commit' });
        return false;
      }
      if (typeof assertCommitted === 'function') {
        const asserted = assertCommitted({ ...sides, direction, prepared });
        if (isPromiseLike(asserted)) {
          throw new Error('Style command post-commit assertion must be synchronous.');
        }
        if (asserted === false) throw new Error('Style command post-commit invariant failed.');
      }
      return true;
    } catch (error) {
      try {
        restore(snapshot, { direction, prepared, reason: 'commit-error', error });
      } catch (rollbackError) {
        const integrityError = new Error(
          `Feature style commit failed and exact rollback also failed: ${rollbackError.message}`
        );
        integrityError.name = 'StyleCommandIntegrityError';
        integrityError.cause = error;
        throw integrityError;
      }
      throw error;
    }
  };

  return Object.freeze({
    label,
    metadata: immutableMetadata,
    counters: immutableCounters,
    apply: ({ direction = 'apply' } = {}) => execute(direction === 'redo' ? 'redo' : 'apply'),
    revert: () => execute('undo'),
    snapshot: ({ direction = 'apply' } = {}) => {
      const snapshot = captureExactState({ direction, compensation: true });
      if (isPromiseLike(snapshot)) {
        throw new Error('Style command compensation capture must be synchronous.');
      }
      return snapshot;
    },
    compensate: ({ direction = 'apply', snapshot, error } = {}) => (
      restore(snapshot, { direction, reason: 'history-compensation', error })
    ),
    estimateBytes: () => {
      // Keep the initially accounted artifact envelope alive for exactly as
      // long as History retains this command.
      void retainedPreparedArtifacts;
      return estimatedBytes;
    }
  });
};
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
