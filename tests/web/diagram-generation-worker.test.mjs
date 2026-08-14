import assert from 'node:assert/strict';

globalThis.self = {};

const {
  collectGenerationResultTransferList,
  resolveGenerationCleanupOutcome,
  serializeError
} = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

const svgBytes = new TextEncoder().encode('<svg />');
const metadataBytes = new TextEncoder().encode('{"featureCatalog":{}}');
assert.deepEqual(
  collectGenerationResultTransferList({
    results: [{ content: svgBytes }, { content: '<svg />' }],
    metadata: metadataBytes
  }),
  [svgBytes.buffer, metadataBytes.buffer]
);

const pythonPayload = {
  error: {
    type: 'ValidationError',
    message: 'primary Python render failure',
    traceback: 'Traceback: primary Python render failure'
  }
};
const pythonResult = resolveGenerationCleanupOutcome({
  result: pythonPayload,
  destroyError: new Error('proxy destroy failed'),
  workspaceError: new Error('workspace invariant failed')
});
assert.equal(pythonResult.error.type, 'ValidationError');
assert.equal(pythonResult.error.message, 'primary Python render failure');
assert.deepEqual(pythonResult.error.notes, [
  'Temporary render handle cleanup also failed: Error: proxy destroy failed',
  'Temporary render workspace cleanup also failed: Error: workspace invariant failed'
]);
assert.match(pythonResult.error.traceback, /proxy destroy failed/);
assert.match(pythonResult.error.traceback, /workspace invariant failed/);

const primaryJsError = new TypeError('primary JavaScript render failure');
assert.throws(
  () => resolveGenerationCleanupOutcome({
    primaryError: primaryJsError,
    destroyError: new Error('proxy destroy failed'),
    workspaceError: new Error('workspace invariant failed')
  }),
  (error) => {
    const serialized = serializeError(error);
    assert.equal(serialized.name, 'TypeError');
    assert.equal(serialized.message, 'TypeError: primary JavaScript render failure');
    assert.match(serialized.stack, /proxy destroy failed/);
    assert.match(serialized.stack, /workspace invariant failed/);
    return true;
  }
);

const destroyError = new Error('proxy destroy failed');
assert.throws(
  () => resolveGenerationCleanupOutcome({
    result: { results: [] },
    destroyError,
    workspaceError: new Error('workspace invariant failed')
  }),
  (error) => {
    assert.equal(error, destroyError);
    assert.match(error.stack, /workspace invariant failed/);
    return true;
  }
);

const workspaceError = new Error('workspace invariant failed');
assert.throws(
  () => resolveGenerationCleanupOutcome({
    result: { results: [] },
    workspaceError
  }),
  (error) => error === workspaceError
);

const structuredError = serializeError({
  type: 'ValidationError',
  message: { summary: 'Annotation target is missing', details: ['row 2'] },
  traceback: 'Traceback (most recent call last): secret input'
});
assert.doesNotMatch(structuredError.message, /\[object Object\]|Traceback/);
assert.match(structuredError.message, /Annotation target is missing/);
