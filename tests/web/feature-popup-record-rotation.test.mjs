import assert from 'node:assert/strict';
import test from 'node:test';

import {
  createFeatureRecordRotationWorkflow
} from '../../gbdraw/web/js/app/record-display/feature-record-rotation.js';
import { resolveFeatureAnchor } from '../../gbdraw/web/js/app/record-display/feature-anchor.js';
import { RECORD_TARGET_NOT_DISCOVERED } from '../../gbdraw/web/js/app/record-display-options.js';

// The reason feature-record-rotation.js gives when the feature's input differs
// from the Result's; the harness action stubs throw it for a stale feature.
const INPUT_DIFFERS = "This feature's input file differs from the one the current Result was drawn from. "
  + 'Generate Diagram to use the current file.';
const feature = {
  record_key: 'record-1',
  biological_feature_id: 'feature-1',
  type: 'CDS'
};
const unrelatedGlobalSelection = [{
  record_key: 'other-record',
  biological_feature_id: 'other-feature'
}];

// The harness resolves through the real resolver: a 1,000 bp feature at
// 2000..2999 on a 9,000 bp circular record.
const createHarness = ({
  strand = '+',
  precision = 'exact',
  committedReverseComplement = false
} = {}) => {
  let applyCount = 0;
  let currentFeature = feature;
  let releaseApply = null;
  let blockApply = false;
  let applyOutcome = { status: 'ok' };
  // The record display draft that Generate would apply, as targetForFeature
  // reports it when it differs from the committed request.
  let pendingTransform = null;
  let stageOutcome = { status: 'ok' };
  const stagedIntents = [];
  const seenFeatures = [];
  const intents = [];
  const action = {
    currentFeature(candidate) {
      if (!currentFeature) throw new Error(INPUT_DIFFERS);
      assert.equal(candidate.record_key, currentFeature.record_key);
      return currentFeature;
    },
    resolve({ feature: explicitFeature, intent }) {
      seenFeatures.push(explicitFeature);
      intents.push(intent);
      if (!currentFeature) throw new Error(INPUT_DIFFERS);
      return {
        feature: currentFeature,
        row: { key: 'row-1' },
        target: {
          recordId: 'NC_000001.1',
          recordKey: 'record-1',
          recordLength: 9000,
          committedReverseComplement,
          pendingTransform
        },
        resolved: resolveFeatureAnchor({
          recordLength: 9000,
          effectiveCircular: true,
          cropped: false,
          currentReverseComplement: committedReverseComplement,
          identity: { recordKey: 'record-1', biologicalFeatureId: 'feature-1' },
          parts: [{ start: 1999, end: 2999, strand: strand === 'unstranded' ? 'undefined' : strand }],
          profile: {
            precision,
            operator: 'single',
            partOrder: strand === 'unstranded' ? 'source-forward' : 'biological',
            strand
          },
          intent
        })
      };
    },
    async apply({ feature: explicitFeature, intent }) {
      applyCount += 1;
      seenFeatures.push(explicitFeature);
      intents.push(intent);
      if (!currentFeature) throw new Error(INPUT_DIFFERS);
      if (blockApply) {
        await new Promise((resolve) => { releaseApply = resolve; });
      }
      if (applyOutcome.status === 'ok') currentFeature = { ...currentFeature, regenerated: true };
      return applyOutcome;
    },
    async stage({ feature: explicitFeature, intent }) {
      seenFeatures.push(explicitFeature);
      stagedIntents.push(intent);
      if (!currentFeature) throw new Error(INPUT_DIFFERS);
      return stageOutcome;
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
    lastIntent: () => intents.at(-1),
    applyCount: () => applyCount,
    stagedIntents,
    setApplyOutcome: (value) => { applyOutcome = value; },
    setStageOutcome: (value) => { stageOutcome = value; },
    setPendingTransform: (value) => { pendingTransform = value; },
    setCurrentFeature: (value) => { currentFeature = value; },
    blockApply: () => { blockApply = true; },
    releaseApply: () => releaseApply?.()
  };
};

const choiceStates = (choices) => choices.map(({ value, enabled, message }) => [value, enabled, message]);

test('popup draft defaults to Start of the record for the explicitly opened feature', () => {
  const harness = createHarness();
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;

  assert.equal(draft.position, 'start');
  assert.equal(draft.reference, 'five-prime');
  assert.equal(draft.offsetText, '0');
  assert.equal(draft.orientForward, false);
  assert.equal(draft.recordLabel, 'NC_000001.1');
  assert.deepEqual(draft.positionChoices.map(({ value, label }) => [value, label]), [
    ['start', 'Start of the record'],
    ['end', 'End of the record'],
    ['custom', 'Custom position']
  ]);
  assert.deepEqual(draft.referenceChoices.map(({ value, label }) => [value, label]), [
    ['five-prime', 'this feature’s 5′ end'],
    ['midpoint', 'this feature’s midpoint'],
    ['three-prime', 'this feature’s 3′ end'],
    ['feature-end', 'just after the feature']
  ]);
  assert.ok(draft.positionChoices.every(({ enabled }) => enabled));
  assert.deepEqual(harness.lastIntent(), {
    placement: 'feature-start', anchor: null, offsetBp: 0, orientForward: false
  });
  assert.equal(draft.startCoordinate, 2000);
  assert.equal(draft.preview, 'NC_000001.1 will start at 2,000 · orientation unchanged');
  assert.equal(draft.canApply, true);
  assert.equal(draft.disabledReason, '');
  assert.equal(harness.seenFeatures.at(-1), feature);
  assert.notEqual(harness.seenFeatures.at(-1), unrelatedGlobalSelection[0]);
});

test('Start of the record puts a minus-strand feature first with or without orientation', () => {
  const harness = createHarness({ strand: '-' });
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  assert.equal(draft.startCoordinate, 2000);
  assert.equal(draft.preview, 'NC_000001.1 will start at 2,000 · orientation unchanged');

  harness.workflow.setOrientForward(true);
  assert.equal(draft.startCoordinate, 2999);
  assert.equal(
    draft.preview,
    'NC_000001.1 will start at 2,999 · reverse-complemented · feature strand − → +'
  );

  const reversed = createHarness({ committedReverseComplement: true });
  reversed.workflow.open({ feature });
  assert.equal(reversed.workflow.draft.startCoordinate, 2999);
  assert.equal(reversed.workflow.draft.preview, 'NC_000001.1 will start at 2,999 · orientation unchanged');
  reversed.workflow.setOrientForward(true);
  assert.equal(reversed.workflow.draft.startCoordinate, 2000);
  assert.equal(
    reversed.workflow.draft.preview,
    'NC_000001.1 will start at 2,000 · no longer reverse-complemented · feature strand − → +'
  );
});

test('End of the record and Custom position map to feature-end and anchor intents', () => {
  const harness = createHarness();
  const { workflow } = harness;
  workflow.open({ feature });

  workflow.setPosition('end');
  assert.deepEqual(harness.lastIntent(), {
    placement: 'feature-end', anchor: null, offsetBp: 0, orientForward: false
  });
  assert.equal(workflow.draft.startCoordinate, 3000);

  workflow.setPosition('custom');
  assert.equal(workflow.draft.reference, 'five-prime');
  for (const [reference, expected] of [
    ['five-prime', { placement: 'anchor', anchor: 'five-prime' }],
    ['midpoint', { placement: 'anchor', anchor: 'midpoint' }],
    ['three-prime', { placement: 'anchor', anchor: 'three-prime' }],
    ['feature-end', { placement: 'feature-end', anchor: null }]
  ]) {
    workflow.setReference(reference);
    workflow.setOffset('-4');
    assert.deepEqual(harness.lastIntent(), { ...expected, offsetBp: -4, orientForward: false });
    assert.equal(workflow.draft.canApply, true, reference);
  }
  assert.equal(workflow.draft.startCoordinate, 2996);

  workflow.setOffset('1.5');
  assert.match(workflow.draft.offsetError, /whole number/);
  assert.equal(workflow.draft.disabledReason, '', 'the offset error is shown once, next to the field');
  assert.equal(workflow.draft.preview, '');
  assert.equal(workflow.draft.canApply, false);
  workflow.setOffset('');
  assert.match(workflow.draft.offsetError, /Enter/);
  workflow.setOffset('9007199254740992');
  assert.match(workflow.draft.offsetError, /supported integer range/);

  workflow.setPosition('start');
  assert.equal(workflow.draft.offsetError, '', 'the hidden offset field adds no reason');
  assert.equal(workflow.draft.canApply, true);
  assert.equal(harness.lastIntent().offsetBp, 0);
});

test('unavailable choices are disabled with their own reason and the selection stays valid', () => {
  const harness = createHarness({ strand: 'unstranded' });
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  assert.equal(draft.targetFailure, '');
  assert.deepEqual(choiceStates(draft.positionChoices), [
    ['start', false, 'Needs a known feature strand.'],
    ['end', true, ''],
    ['custom', true, '']
  ]);
  assert.equal(draft.position, 'end', 'the first available choice replaces an unavailable one');
  assert.equal(draft.canApply, true);
  assert.equal(draft.orientCapability.enabled, false);
  assert.match(draft.orientCapability.message, /known feature strand/);
  assert.deepEqual(draft.referenceChoices.map(({ value, enabled }) => [value, enabled]), [
    ['five-prime', false], ['midpoint', true], ['three-prime', false], ['feature-end', true]
  ]);
  assert.match(draft.referenceChoices[0].message, /known strand/);

  harness.workflow.setPosition('start');
  assert.equal(draft.position, 'end');
  harness.workflow.setPosition('custom');
  assert.equal(draft.reference, 'midpoint');
  assert.equal(harness.lastIntent().anchor, 'midpoint');
  assert.equal(draft.canApply, true);
});

test('a location limit that disables every placement is one whole-target reason', () => {
  const harness = createHarness({ precision: 'fuzzy' });
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  assert.equal(draft.targetFailure, 'location-fuzzy');
  assert.equal(draft.disabledReason, 'This operation requires exact feature coordinates.');
  assert.ok([...draft.positionChoices, ...draft.referenceChoices]
    .every(({ enabled, message }) => !enabled && message === ''));
  assert.equal(draft.orientCapability.message, '');
  assert.equal(draft.canApply, false);
  assert.equal(draft.preview, '');
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
  assert.match(harness.workflow.draft.disabledReason, /input file differs/);
  // One reason source: the reason line states it; the status line does not repeat it.
  assert.equal(harness.workflow.draft.status, '');
  assert.equal(harness.applyCount(), 1);
});

test('a failed Apply and regenerate states its reason once', async () => {
  const harness = createHarness();
  harness.workflow.open({ feature });
  harness.action.apply = async () => {
    harness.setCurrentFeature(null);
    throw new Error(INPUT_DIFFERS);
  };
  const outcome = await harness.workflow.apply();
  const { draft } = harness.workflow;
  assert.equal(outcome.status, 'error');
  assert.equal(draft.disabledReason, INPUT_DIFFERS);
  assert.equal(draft.status, 'Record rotation failed. The previous Result was kept.');
  assert.equal(draft.statusKind, 'error');

  const transient = createHarness();
  transient.workflow.open({ feature });
  transient.action.apply = async () => { throw new Error('Worker stopped.'); };
  await transient.workflow.apply();
  assert.equal(transient.workflow.draft.disabledReason, '');
  assert.equal(transient.workflow.draft.status, 'Worker stopped. The previous Result was kept.');
});

test('Apply on Generate stages the resolved transform with the same validation and no candidate', async () => {
  const harness = createHarness();
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  harness.workflow.setPosition('end');

  const outcome = await harness.workflow.stage();
  assert.equal(outcome.status, 'ok');
  assert.equal(harness.applyCount(), 0, 'staging runs no candidate');
  assert.deepEqual(harness.stagedIntents, [{
    placement: 'feature-end', anchor: null, offsetBp: 0, orientForward: false
  }]);
  assert.equal(harness.seenFeatures.at(-1), feature);
  assert.equal(draft.status, 'Record rotation will apply on the next Generate Diagram.');
  assert.equal(draft.statusKind, 'success');
  assert.equal(draft.pending, false);
  assert.equal(draft.canApply, true);

  // Both apply buttons share canApply: an invalid offset stages nothing.
  harness.workflow.setPosition('custom');
  harness.workflow.setOffset('1.5');
  assert.equal(draft.canApply, false);
  assert.deepEqual(await harness.workflow.stage(), { status: 'disabled' });
  assert.equal(harness.stagedIntents.length, 1);
  assert.equal(draft.status, '', 'the offset error is shown once, next to the field');

  // A target that went stale before the click: one reason, no stage.
  harness.workflow.setOffset('0');
  harness.setCurrentFeature(null);
  assert.deepEqual(await harness.workflow.stage(), { status: 'disabled' });
  assert.match(draft.disabledReason, /input file differs/);
  assert.equal(draft.status, '');
  assert.equal(harness.stagedIntents.length, 1);
});

test('Apply on Generate reports a busy Session and is disabled while regenerating', async () => {
  const harness = createHarness();
  harness.workflow.open({ feature });
  harness.setStageOutcome({ status: 'busy', reason: 'Saving session. Retry after saving finishes.' });
  assert.equal((await harness.workflow.stage()).status, 'busy');
  assert.equal(harness.workflow.draft.status, 'Saving session. Retry after saving finishes.');
  assert.equal(harness.workflow.draft.statusKind, 'error');

  harness.setStageOutcome({ status: 'ok' });
  harness.blockApply();
  const regenerating = harness.workflow.apply();
  assert.equal(harness.workflow.draft.canApply, false);
  assert.deepEqual(await harness.workflow.stage(), { status: 'pending' });
  assert.equal(harness.stagedIntents.length, 1);
  harness.releaseApply();
  await regenerating;
});

test('a pending draft for the record is shown when the popup opens', () => {
  const harness = createHarness();
  harness.setPendingTransform({ startCoordinate: 4500, reverseComplement: true });
  harness.workflow.open({ feature });
  const { draft } = harness.workflow;
  assert.equal(
    draft.pendingPreview,
    'Pending for Generate: NC_000001.1 will start at 4,500 · reverse-complemented'
  );
  assert.equal(draft.preview, 'NC_000001.1 will start at 2,000 · orientation unchanged');

  // A pending draft equal to the chosen outcome is stated once.
  harness.setPendingTransform({ startCoordinate: 2000, reverseComplement: false });
  harness.workflow.recompute();
  assert.equal(
    draft.pendingPreview,
    'Pending for Generate: NC_000001.1 will start at 2,000 · orientation unchanged'
  );
  assert.equal(draft.preview, '');
  assert.equal(draft.canApply, true);

  // An orientation-only draft starts at the record's default start.
  harness.setPendingTransform({ startCoordinate: null, reverseComplement: true });
  harness.workflow.recompute();
  assert.equal(
    draft.pendingPreview,
    'Pending for Generate: NC_000001.1 will start at 9,000 · reverse-complemented'
  );

  harness.setPendingTransform(null);
  harness.workflow.recompute();
  assert.equal(draft.pendingPreview, '');
  assert.equal(draft.preview, 'NC_000001.1 will start at 2,000 · orientation unchanged');
});

for (const change of ['close', 'retarget']) {
  test(`an Apply and regenerate that settles after a popup ${change} writes nothing into the current draft (F7)`, async () => {
    const harness = createHarness();
    harness.workflow.open({ feature });
    harness.blockApply();
    const regenerating = harness.workflow.apply();
    harness.workflow.close();
    if (change === 'retarget') {
      harness.workflow.open({ feature: { ...feature, biological_feature_id: 'feature-2' } });
    }
    const current = JSON.parse(JSON.stringify(harness.workflow.draft));
    harness.releaseApply();
    assert.equal((await regenerating).status, 'ok');
    assert.deepEqual(JSON.parse(JSON.stringify(harness.workflow.draft)), current);
    assert.deepEqual(harness.rebound, []);
  });
}

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
const enabled = { enabled: true, code: null, message: '' };
const availableCapabilities = {
  featureStart: enabled,
  anchors: { 'five-prime': enabled, midpoint: enabled, 'three-prime': enabled },
  offset: enabled,
  orientForward: enabled,
  featureEnd: enabled
};

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
          capabilities: availableCapabilities,
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
  assert.ok([...draft.positionChoices, ...draft.referenceChoices]
    .every(({ message }) => message === ''));
  assert.equal(draft.orientCapability.message, '');
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
            featureStart: limit,
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
  workflow.setPosition('custom');
  workflow.setOffset('x');
  assert.equal(draft.disabledReason, limit.message, 'the hidden offset field adds no reason');
  assert.equal(draft.offsetError, '');
  assert.deepEqual(choiceStates(draft.positionChoices),
    [['start', false, ''], ['end', false, ''], ['custom', false, '']]);
  assert.ok(draft.referenceChoices.every(({ enabled, message }) => !enabled && message === ''));
  assert.equal(draft.orientCapability.message, '');
  assert.equal(draft.canApply, false);
});
