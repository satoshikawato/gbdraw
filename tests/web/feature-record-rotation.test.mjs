import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import { createFeatureRecordRotationAction } from '../../gbdraw/web/js/app/record-display/feature-record-rotation.js';
import { projectCommittedRecordTransform } from '../../gbdraw/web/js/services/session-request.js';
import { validateAnchorIntent } from '../../gbdraw/web/js/app/record-display-options.js';

const committed = JSON.parse(await readFile(
  'docs/images/h-cli-12/cli_session.json',
  'utf8'
));
const canonicalRecord = committed.renderRequest.records[0];
const catalogFeatures = committed.editorState.featureCatalog.items[0].biologicalFeatures;
const sourceFeature = catalogFeatures.find((feature) => (
  feature.anchorProfile?.precision === 'exact'
  && feature.anchorProfile?.strand === '+'
  && Number.isSafeInteger(feature.start)
  && Number.isSafeInteger(feature.end)
));
const feature = {
  record_key: canonicalRecord.recordKey,
  biological_feature_id: sourceFeature.biologicalFeatureId,
  start: sourceFeature.start,
  end: sourceFeature.end,
  strand: sourceFeature.strand,
  location_parts: [{
    start: sourceFeature.start,
    end: sourceFeature.end,
    strand: sourceFeature.strand,
    display: `${sourceFeature.start + 1}..${sourceFeature.end}`
  }],
  anchorProfile: sourceFeature.anchorProfile
};
const row = { key: 'target-row' };
const target = {
  scope: committed.renderRequest.mode,
  sourceUid: 'circular',
  selector: '#1',
  recordId: sourceFeature.record_id,
  recordLength: 16569,
  recordKey: canonicalRecord.recordKey,
  canonicalRecordKey: canonicalRecord.recordKey,
  source: structuredClone(canonicalRecord.source),
  effectiveCircular: true,
  cropped: false,
  committedReverseComplement: false,
  members: [{
    canonicalRecordKey: canonicalRecord.recordKey,
    recordKey: canonicalRecord.recordKey,
    selector: '#1',
    recordId: sourceFeature.record_id,
    recordLength: 16569,
    committedDisplay: structuredClone(canonicalRecord.display),
    committedReverseComplement: false
  }]
};

test('explicit popup feature coordinates resolver, projector, runner, and target draft owner', async () => {
  let targetLookups = 0;
  let committedTransform = null;
  let executedCanonical = null;
  const beforeCheckpoint = { key: 'target-row', index: -1, draft: null };
  const controls = {
    targetForFeature(clickedFeature) {
      targetLookups += 1;
      assert.equal(clickedFeature, feature);
      return { row, target };
    },
    captureTargetDraft(clickedRow) {
      assert.equal(clickedRow, row);
      return beforeCheckpoint;
    },
    restoreTargetDraft() {
      throw new Error('successful action must not restore the before draft');
    },
    commitResolvedTransform(clickedRow, transform) {
      assert.equal(clickedRow, row);
      committedTransform = structuredClone(transform);
    }
  };
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: controls,
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async (options) => {
      executedCanonical = options.canonical;
      assert.equal(options.captureIntentCheckpoint(), beforeCheckpoint);
      await options.commitIntent();
      return { status: 'ok' };
    }
  });
  const outcome = await action.apply({
    feature,
    intent: {
      placement: 'anchor',
      anchor: 'five-prime',
      offsetBp: -100,
      orientForward: true
    }
  });
  assert.equal(targetLookups, 1);
  assert.equal(outcome.status, 'ok');
  assert.equal(outcome.receipt.recordKey, canonicalRecord.recordKey);
  assert.equal(
    executedCanonical.renderRequest.records[0].display.startCoordinate,
    outcome.resolved.startCoordinate
  );
  assert.equal(
    executedCanonical.renderRequest.records[0].presentation.reverseComplement,
    outcome.resolved.reverseComplement
  );
  assert.deepEqual(committedTransform.anchorIntent, outcome.resolved.provenance);
  assert.equal(committedTransform.startCoordinate, outcome.resolved.startCoordinate);
  assert.equal(committedTransform.reverseComplement, outcome.resolved.reverseComplement);
});

test('Start of the record commits the resolved anchor in the persisted provenance format', async () => {
  let committedTransform = null;
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: {
      targetForFeature: () => ({ row, target }),
      captureTargetDraft: () => null,
      restoreTargetDraft() { throw new Error('unreachable'); },
      commitResolvedTransform(clickedRow, transform) {
        committedTransform = structuredClone(transform);
      }
    },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async (options) => {
      await options.commitIntent();
      return { status: 'ok' };
    }
  });
  const outcome = await action.apply({
    feature,
    intent: { placement: 'feature-start', anchor: null, offsetBp: 0, orientForward: false }
  });
  assert.equal(outcome.status, 'ok');
  assert.equal(committedTransform.startCoordinate, sourceFeature.start + 1);
  assert.equal(committedTransform.reverseComplement, false);
  assert.equal(validateAnchorIntent(committedTransform.anchorIntent), committedTransform.anchorIntent);
  assert.deepEqual(committedTransform.anchorIntent, {
    schema: 1,
    recordKey: canonicalRecord.recordKey,
    biologicalFeatureId: sourceFeature.biologicalFeatureId,
    placement: 'anchor',
    anchor: 'five-prime',
    offsetBp: 0,
    orientForward: false
  });
});

test('disabled source profile rejects before candidate execution', async () => {
  let executions = 0;
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: { targetForFeature: () => ({ row, target }) },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async () => { executions += 1; }
  });
  await assert.rejects(action.apply({
    feature: {
      ...feature,
      anchorProfile: {
        precision: 'unavailable',
        operator: 'unknown',
        partOrder: 'ambiguous',
        strand: '+'
      }
    },
    intent: {
      placement: 'anchor', anchor: 'five-prime', offsetBp: 0, orientForward: false
    }
  }), /Generate again/);
  assert.equal(executions, 0);
});

test('Apply re-resolves the stable popup identity and rejects a stale feature', async () => {
  let currentFeature = feature;
  let executions = 0;
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: { targetForFeature: () => ({ row, target }) },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    resolveCurrentFeature: () => currentFeature,
    isCurrentFeature: () => true,
    runCommittedCanonicalCandidate: async () => {
      executions += 1;
      return { status: 'ok' };
    }
  });
  assert.equal(action.resolve({
    feature,
    intent: {
      placement: 'anchor', anchor: 'five-prime', offsetBp: 0, orientForward: false
    }
  }).resolved.eligibility.enabled, true);

  currentFeature = null;
  await assert.rejects(action.apply({
    feature,
    intent: {
      placement: 'anchor', anchor: 'five-prime', offsetBp: 0, orientForward: false
    }
  }), /no longer present/);
  assert.equal(executions, 0);
});

test('a feature whose input differs from the Result input states the reason and the next step', async () => {
  let executions = 0;
  let staged = 0;
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: {
      targetForFeature: () => ({ row, target }),
      setResolvedTransform: async () => { staged += 1; }
    },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    resolveCurrentFeature: () => feature,
    isCurrentFeature: () => false,
    runCommittedCanonicalCandidate: async () => {
      executions += 1;
      return { status: 'ok' };
    }
  });
  const reason = "This feature's input file differs from the one the current Result was drawn from. "
    + 'Generate Diagram to use the current file.';
  const intent = { placement: 'anchor', anchor: 'five-prime', offsetBp: 0, orientForward: false };
  assert.throws(() => action.resolve({ feature, intent }), { message: reason });
  await assert.rejects(action.apply({ feature, intent }), { message: reason });
  await assert.rejects(action.stage({ feature, intent }), { message: reason });
  assert.equal(executions, 0);
  assert.equal(staged, 0);
});

test('the record read uses the injected mode discovery and has no fallback reader', async () => {
  let reads = 0;
  const owners = {
    recordDisplayControls: { targetForFeature: () => ({ row, target }) },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async () => ({ status: 'ok' })
  };
  const action = createFeatureRecordRotationAction({
    ...owners,
    readRecords: async () => { reads += 1; return { status: 'ok' }; }
  });
  assert.deepEqual(await action.readRecords(), { status: 'ok' });
  assert.equal(reads, 1);
  assert.equal((await createFeatureRecordRotationAction(owners).readRecords()).status, 'unavailable');
});

test('Apply on Generate writes the resolved transform as one History step and runs no candidate', async () => {
  const staged = [];
  let executions = 0;
  const action = createFeatureRecordRotationAction({
    recordDisplayControls: {
      targetForFeature: () => ({ row, target }),
      setResolvedTransform: async (clickedRow, transform) => {
        assert.equal(clickedRow, row);
        staged.push(structuredClone(transform));
      },
      commitResolvedTransform() { throw new Error('staging uses the undoable writer'); }
    },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async () => { executions += 1; return { status: 'ok' }; }
  });
  const outcome = await action.stage({
    feature,
    intent: { placement: 'feature-end', anchor: null, offsetBp: 0, orientForward: false }
  });
  assert.equal(outcome.status, 'ok');
  assert.equal(executions, 0);
  assert.deepEqual(staged, [{
    startCoordinate: sourceFeature.end + 1,
    reverseComplement: false,
    anchorIntent: outcome.resolved.provenance
  }]);
  assert.equal(outcome.resolved.provenance.placement, 'feature-end');

  const busy = { status: 'busy', reason: 'Saving session. Retry after saving finishes.' };
  const busyAction = createFeatureRecordRotationAction({
    recordDisplayControls: {
      targetForFeature: () => ({ row, target }),
      setResolvedTransform: async () => busy
    },
    getCommittedSession: () => committed,
    projectCommittedRecordTransform,
    runCommittedCanonicalCandidate: async () => { executions += 1; }
  });
  assert.equal(await busyAction.stage({
    feature,
    intent: { placement: 'feature-start', anchor: null, offsetBp: 0, orientForward: false }
  }), busy);
  await assert.rejects(busyAction.stage({
    feature: { ...feature, anchorProfile: { ...feature.anchorProfile, precision: 'unavailable' } },
    intent: { placement: 'feature-start', anchor: null, offsetBp: 0, orientForward: false }
  }), /Generate again/);
  assert.equal(executions, 0);
});
