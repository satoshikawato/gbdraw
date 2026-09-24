import assert from 'node:assert/strict';
import test from 'node:test';

import {
  createFeatureRecordRotationWorkflow
} from '../../gbdraw/web/js/app/record-display/feature-record-rotation.js';

const feature = {
  record_key: 'record-1',
  biological_feature_id: 'feature-1',
  type: 'CDS'
};
const unrelatedGlobalSelection = [{
  record_key: 'other-record',
  biological_feature_id: 'other-feature'
}];

const capabilities = ({ stranded = true, available = true } = {}) => ({
  anchors: {
    'five-prime': {
      enabled: available && stranded,
      message: stranded ? '' : 'Biological ends require a known strand.'
    },
    midpoint: { enabled: available, message: available ? '' : 'Location is ambiguous.' },
    'three-prime': {
      enabled: available && stranded,
      message: stranded ? '' : 'Biological ends require a known strand.'
    }
  },
  offset: { enabled: available, message: available ? '' : 'Location is ambiguous.' },
  orientForward: {
    enabled: available && stranded,
    message: stranded ? '' : 'Orient forward requires a known strand.'
  },
  featureEnd: { enabled: available, message: available ? '' : 'Traversal is ambiguous.' }
});

const createHarness = ({ stranded = true, available = true } = {}) => {
  let applyCount = 0;
  let currentFeature = feature;
  let releaseApply = null;
  let blockApply = false;
  const seenFeatures = [];
  const action = {
    currentFeature(candidate) {
      if (!currentFeature) throw new Error('The popup feature source changed.');
      assert.equal(candidate.record_key, currentFeature.record_key);
      return currentFeature;
    },
    resolve({ feature: explicitFeature, intent }) {
      seenFeatures.push(explicitFeature);
      if (!currentFeature) throw new Error('The popup feature source changed.');
      const offsetValid = Number.isSafeInteger(intent.offsetBp);
      const allowed = available && offsetValid
        && (intent.placement === 'feature-end'
          || intent.anchor === 'midpoint'
          || stranded);
      return {
        feature: currentFeature,
        row: { key: 'row-1' },
        target: {
          recordId: 'NC_000001.1',
          recordKey: 'record-1',
          committedReverseComplement: false
        },
        resolved: {
          capabilities: capabilities({ stranded, available }),
          eligibility: {
            enabled: allowed,
            message: allowed ? ''
              : !offsetValid ? 'Offset must be an integer.'
                : stranded ? 'Location is ambiguous.' : 'Biological ends require a known strand.'
          },
          startCoordinate: allowed ? 11 + intent.offsetBp : null,
          reverseComplement: allowed ? Boolean(intent.orientForward && !stranded) : null,
          displayedStrand: stranded ? { before: '+', after: '+' } : {
            before: 'unstranded', after: 'unstranded'
          }
        }
      };
    },
    async apply({ feature: explicitFeature }) {
      applyCount += 1;
      seenFeatures.push(explicitFeature);
      if (!currentFeature) throw new Error('The popup feature source changed.');
      if (blockApply) {
        await new Promise((resolve) => { releaseApply = resolve; });
      }
      currentFeature = { ...currentFeature, regenerated: true };
      return { status: 'ok' };
    }
  };
  const rebound = [];
  const workflow = createFeatureRecordRotationWorkflow({
    action,
    onRebind: (next) => rebound.push(next)
  });
  return {
    action,
    workflow,
    rebound,
    seenFeatures,
    applyCount: () => applyCount,
    setCurrentFeature: (value) => { currentFeature = value; },
    blockApply: () => { blockApply = true; },
    releaseApply: () => releaseApply?.()
  };
};

test('popup draft defaults and previews use only the explicitly opened feature', () => {
  const harness = createHarness();
  harness.workflow.open({ feature, featureLabel: 'target CDS' });

  assert.equal(harness.workflow.draft.anchor, 'five-prime');
  assert.equal(harness.workflow.draft.offsetText, '0');
  assert.equal(harness.workflow.draft.orientForward, false);
  assert.equal(harness.workflow.draft.recordLabel, 'NC_000001.1');
  assert.equal(harness.workflow.draft.featureLabel, 'target CDS');
  assert.equal(harness.workflow.draft.startCoordinate, 11);
  assert.equal(harness.seenFeatures.at(-1), feature);
  assert.notEqual(harness.seenFeatures.at(-1), unrelatedGlobalSelection[0]);

  harness.workflow.setAnchor('midpoint');
  harness.workflow.setOffset('-4');
  assert.equal(harness.workflow.draft.startCoordinate, 7);
  assert.equal(harness.workflow.draft.canApply, true);

  harness.workflow.setOffset('1.5');
  assert.match(harness.workflow.draft.offsetError, /whole number/);
  assert.equal(harness.workflow.draft.canApply, false);
  harness.workflow.setOffset('');
  assert.match(harness.workflow.draft.offsetError, /Enter/);
  harness.workflow.setOffset('9007199254740992');
  assert.match(harness.workflow.draft.offsetError, /supported integer range/);
});

test('feature-end preset resets offset and later edits become custom placement', () => {
  const { workflow } = createHarness();
  workflow.open({ feature });
  workflow.setOffset('12');
  workflow.placeAtFeatureEnd();
  assert.equal(workflow.draft.placement, 'feature-end');
  assert.equal(workflow.draft.offsetText, '0');
  assert.equal(workflow.draft.exactFeatureEnd, true);
  assert.equal(workflow.draft.placementLabel, 'Exactly at feature end');

  workflow.setOffset('3');
  assert.equal(workflow.draft.exactFeatureEnd, false);
  assert.equal(workflow.draft.placementLabel, 'Feature end with custom offset');
});

test('operation capabilities keep safe unstranded midpoint and feature-end choices', () => {
  const { workflow } = createHarness({ stranded: false });
  workflow.open({ feature });
  assert.equal(workflow.draft.canApply, false);
  assert.match(workflow.draft.disabledReason, /known strand/);
  assert.equal(
    workflow.draft.anchorChoices.find(({ value }) => value === 'midpoint').enabled,
    true
  );
  assert.equal(workflow.draft.orientCapability.enabled, false);
  assert.match(workflow.draft.orientCapability.message, /known strand/);

  workflow.setAnchor('midpoint');
  assert.equal(workflow.draft.canApply, true);
  workflow.placeAtFeatureEnd();
  assert.equal(workflow.draft.canApply, true);
});

test('Apply rechecks freshness, suppresses double submit, and rebinds stable identity', async () => {
  const harness = createHarness();
  harness.workflow.open({ feature });
  harness.blockApply();
  const first = harness.workflow.apply();
  const second = await harness.workflow.apply();
  assert.deepEqual(second, { status: 'pending' });
  assert.equal(harness.applyCount(), 1);
  harness.releaseApply();
  assert.equal((await first).status, 'ok');
  assert.equal(harness.rebound.length, 1);
  assert.equal(harness.rebound[0].regenerated, true);
  assert.match(harness.workflow.draft.status, /applied and regenerated/);

  harness.setCurrentFeature(null);
  const stale = await harness.workflow.apply();
  assert.equal(stale.status, 'disabled');
  assert.match(harness.workflow.draft.disabledReason, /source changed/);
  assert.match(harness.workflow.draft.status, /source changed/);
  assert.equal(harness.workflow.draft.statusKind, 'error');
  assert.equal(harness.applyCount(), 1);
});

test('Cancel discards only the ephemeral popup draft', () => {
  const harness = createHarness();
  const committedSidebarState = { startCoordinate: 1, reverseComplement: false };
  const currentResult = { name: 'before.svg' };
  harness.workflow.open({ feature });
  harness.workflow.setOffset('25');
  harness.workflow.cancel();

  assert.equal(harness.workflow.draft.active, false);
  assert.equal(harness.workflow.draft.feature, null);
  assert.equal(harness.applyCount(), 0);
  assert.deepEqual(committedSidebarState, { startCoordinate: 1, reverseComplement: false });
  assert.deepEqual(currentResult, { name: 'before.svg' });
});
