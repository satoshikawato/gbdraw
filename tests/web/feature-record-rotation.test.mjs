import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import { createFeatureRecordRotationAction } from '../../gbdraw/web/js/app/record-display/feature-record-rotation.js';
import { projectCommittedRecordTransform } from '../../gbdraw/web/js/services/session-request.js';

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
