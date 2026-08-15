import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';

globalThis.self = {};
if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  buildGeneratedArtifactTransportIdentity,
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
const firstIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes
});
const sameIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: metadataBytes
});
const changedIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'changed.svg', content: svgBytes }],
  metadata: metadataBytes
});
const changedContentIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: new TextEncoder().encode('<svg>b</svg>') }],
  metadata: metadataBytes
});
const changedMetadataIdentity = await buildGeneratedArtifactTransportIdentity({
  results: [{ name: 'out.svg', content: svgBytes }],
  metadata: new TextEncoder().encode('{"featureCatalog":{"schema":4}}')
});
assert.match(firstIdentity.fingerprint, /^[0-9a-f]{64}$/);
assert.equal(firstIdentity.fingerprint, sameIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedContentIdentity.fingerprint);
assert.notEqual(firstIdentity.fingerprint, changedMetadataIdentity.fingerprint);
assert.equal(
  firstIdentity.retainedBytes,
  svgBytes.byteLength
    + metadataBytes.byteLength * 2
    + new TextEncoder().encode('out.svg').byteLength
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
