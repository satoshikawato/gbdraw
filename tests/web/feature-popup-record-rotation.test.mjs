import assert from 'node:assert/strict';
import test from 'node:test';

import {
  createFeatureRecordRotationWorkflow
} from '../../gbdraw/web/js/app/record-display/feature-record-rotation.js';
import { RECORD_TARGET_NOT_DISCOVERED } from '../../gbdraw/web/js/app/record-display-options.js';

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

const targetFailure = (kind, message) => Object.assign(new Error(message), { kind });

const createDiscoveryHarness = ({ readResult = { discovered: true } } = {}) => {
  let discovered = false;
  let readError = '';
  let reads = 0;
  let releaseRead = null;
  const action = {
    currentFeature: (candidate) => candidate,
    resolve() {
      if (readError) throw targetFailure('discovery-failed', readError);
      if (!discovered) {
        throw targetFailure(RECORD_TARGET_NOT_DISCOVERED, 'Records for this feature are not read yet.');
      }
      return {
        feature,
        row: { key: 'row-1' },
        target: { recordId: 'NC_000001.1', recordKey: 'record-1', committedReverseComplement: false },
        resolved: {
          capabilities: capabilities(),
          eligibility: { enabled: true, message: '' },
          startCoordinate: 11,
          reverseComplement: false,
          displayedStrand: { before: '+', after: '+' }
        }
      };
    },
    async apply() { throw new Error('Apply is not part of record reading.'); },
    async readRecords() {
      reads += 1;
      await new Promise((resolve) => { releaseRead = resolve; });
      if (readResult.discovered) discovered = true;
      if (readResult.error) readError = readResult.error;
      return readResult.outcome;
    }
  };
  const workflow = createFeatureRecordRotationWorkflow({ action });
  return { workflow, reads: () => reads, releaseRead: () => releaseRead?.() };
};

test('unread records after Session Load are read once, then the draft is recomputed', async () => {
  const harness = createDiscoveryHarness();
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  assert.equal(draft.targetFailure, RECORD_TARGET_NOT_DISCOVERED);
  assert.doesNotMatch(draft.disabledReason, /stale|ambiguous/);
  // A whole-target failure is one reason, not one copy per control (F4).
  assert.deepEqual(draft.anchorChoices.map(({ message }) => message), ['', '', '']);
  assert.equal(draft.orientCapability.message, '');
  assert.equal(draft.featureEndCapability.message, '');
  assert.equal(draft.canApply, false);

  const reading = harness.workflow.readRecords();
  assert.equal(draft.reading, true);
  assert.equal(draft.status, 'Reading records…');
  assert.equal(draft.disabledReason, '');
  assert.equal(draft.canApply, false);
  harness.workflow.readRecords();
  assert.equal(harness.reads(), 1);

  harness.releaseRead();
  await reading;
  assert.equal(draft.reading, false);
  assert.equal(draft.targetFailure, '');
  assert.equal(draft.status, '');
  assert.equal(draft.disabledReason, '');
  assert.equal(draft.recordLabel, 'NC_000001.1');
  assert.equal(draft.startCoordinate, 11);
  assert.equal(draft.canApply, true);
  await harness.workflow.readRecords();
  assert.equal(harness.reads(), 1);
});

test('a failed record read shows its own reason and is not retried automatically', async () => {
  const harness = createDiscoveryHarness({
    readResult: { discovered: false, error: 'Records could not be loaded: malformed LOCUS line.' }
  });
  harness.workflow.open({ feature });
  const reading = harness.workflow.readRecords();
  harness.releaseRead();
  await reading;
  const { draft } = harness.workflow;
  assert.equal(draft.targetFailure, 'discovery-failed');
  assert.equal(draft.disabledReason, 'Records could not be loaded: malformed LOCUS line.');
  assert.equal(draft.status, '');
  await harness.workflow.readRecords();
  assert.equal(harness.reads(), 1);
});

test('a busy Session operation leaves records unread with its reason', async () => {
  const harness = createDiscoveryHarness({
    readResult: { discovered: false, outcome: { status: 'busy', reason: 'Saving session. Retry after saving finishes.' } }
  });
  harness.workflow.open({ feature });
  const reading = harness.workflow.readRecords();
  harness.releaseRead();
  await reading;
  assert.equal(harness.workflow.draft.targetFailure, RECORD_TARGET_NOT_DISCOVERED);
  assert.equal(harness.workflow.draft.disabledReason, 'Saving session. Retry after saving finishes.');
});

test('a record read that settles after the popup closed writes nothing into the next draft', async () => {
  const harness = createDiscoveryHarness();
  harness.workflow.open({ feature });
  const reading = harness.workflow.readRecords();
  harness.workflow.close();
  harness.workflow.open({ feature: { ...feature, biological_feature_id: 'feature-2' } });
  const next = { ...harness.workflow.draft };
  harness.releaseRead();
  await reading;
  assert.equal(harness.workflow.draft.status, next.status);
  assert.equal(harness.workflow.draft.reading, false);
  assert.equal(harness.workflow.draft.targetFailure, RECORD_TARGET_NOT_DISCOVERED);
  assert.equal(harness.workflow.draft.identity.biologicalFeatureId, 'feature-2');
});

test('a limit that disables every control is one whole-target reason', () => {
  const limit = {
    enabled: false,
    code: 'record-not-circular',
    message: 'Record rotation requires an effectively circular record.'
  };
  const workflow = createFeatureRecordRotationWorkflow({
    action: {
      currentFeature: (candidate) => candidate,
      resolve: () => ({
        feature,
        row: { key: 'row-1' },
        target: { recordId: 'NC_001416.1', recordKey: 'record-1', committedReverseComplement: false },
        resolved: {
          capabilities: {
            anchors: { 'five-prime': limit, midpoint: limit, 'three-prime': limit },
            offset: limit,
            orientForward: limit,
            featureEnd: limit
          },
          eligibility: limit,
          startCoordinate: null,
          reverseComplement: null,
          displayedStrand: { before: null, after: null }
        }
      }),
      async apply() { throw new Error('unreachable'); }
    }
  });
  workflow.open({ feature });
  const { draft } = workflow;
  assert.equal(draft.targetFailure, 'record-not-circular');
  assert.equal(draft.recordLabel, 'NC_001416.1');
  assert.equal(draft.disabledReason, limit.message);
  workflow.setOffset('x');
  assert.equal(draft.disabledReason, limit.message, 'the hidden offset field adds no reason');
  assert.deepEqual(draft.anchorChoices.map(({ enabled, message }) => [enabled, message]),
    [[false, ''], [false, ''], [false, '']]);
  assert.deepEqual([draft.orientCapability.message, draft.featureEndCapability.message], ['', '']);
  assert.equal(draft.canApply, false);
});
