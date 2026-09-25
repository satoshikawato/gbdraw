import assert from 'node:assert/strict';
import test from 'node:test';

const ref = (value) => ({ value });
globalThis.window = {
  Vue: { ref, computed: (getter) => ({ get value() { return getter(); } }) }
};
const { createSimilarityAlignmentActions } = await import(
  '../../gbdraw/web/js/app/similarity-alignment.js'
);

const anchor = (recordKey, biologicalFeatureId, sourceFeatureIndex) => ({
  recordKey, biologicalFeatureId, sourceFeatureIndex,
  stableFeatureSvgId: `stable-${biologicalFeatureId}`
});
const member = (recordKey, id, index, strand, representative = false) => ({
  ...anchor(recordKey, id, index), recordIndex: recordKey === 'a' ? 0 : recordKey === 'b' ? 1 : 2,
  featureIndex: index, start: 10 + index * 10, end: 30 + index * 10,
  strand, proteinId: id, sourceProteinId: id, representative, role: 'member',
  gene: `gene-${id}`
});
const reference = member('a', 'clicked', 4, '+');
const otherReference = member('a', 'other', 7, '+', true);
const targetA = member('b', 'target-a', 1, '-', true);
const targetB = member('b', 'target-b', 2, '+');
const targetC = member('c', 'target-c', 3, null);
const renderRequest = {
  mode: 'linear', records: ['a', 'b', 'c'].map((recordKey) => ({
    recordKey, region: null, presentation: { reverseComplement: false }
  }))
};
const state = () => ({
  featureCatalog: ref({ items: [{
    recordKeys: ['a', 'b', 'c'],
    sequenceSources: [{ sequence: 'A'.repeat(300) }, { sequence: 'C'.repeat(300) },
      { sequence: 'G'.repeat(300) }]
  }] }),
  selectedOrthogroupId: ref('og-1'),
  similarityAlignmentPlan: ref(null),
  linearRecordTranslations: ref([]),
  linearSeqs: ['a', 'b', 'c'].map((uid) => ({ uid, region_reverse: false })),
  results: ref([{ name: 'prior.svg' }])
});
const deferred = () => {
  let resolve;
  const promise = new Promise((done) => { resolve = done; });
  return { promise, resolve };
};
const serializedCandidate = (item, request) => {
  const displayedStrand = item.sourceStrand;
  const referenceStrand = request.members.find(({ anchor: selected }) => (
    selected.recordKey === request.reference.recordKey
    && selected.biologicalFeatureId === request.reference.biologicalFeatureId
  ))?.sourceStrand;
  const opposite = displayedStrand !== null && referenceStrand !== null
    && displayedStrand !== referenceStrand;
  return {
    anchor: item.anchor, displayName: item.anchor.biologicalFeatureId,
    sourceStart: item.sourceStart, sourceEnd: item.sourceEnd,
    displayCenter: (item.sourceStart + item.sourceEnd) / 2,
    displayedStrand, hidden: false, representative: item.representative,
    role: item.role, usable: true, directEvidence: [],
    orientation: {
      preserve: { effect: 'preserve', effectiveReverseComplement: false },
      match_reference: {
        effect: displayedStrand === null ? 'preserve_unknown_strand'
          : opposite ? 'reverse_whole_record' : 'preserve',
        effectiveReverseComplement: opposite
      }
    }
  };
};
const responseFor = (request, { ambiguous = false, missing = false } = {}) => {
  const records = request.records.map(({ recordKey }) => {
    if (recordKey === request.reference.recordKey) return {
      kind: 'decision', recordKey, status: 'reference', rationale: 'reference',
      reviewReason: null, anchor: request.reference,
      orientationPolicy: 'preserve', effectiveReverseComplement: null, candidates: []
    };
    const candidates = request.members.filter((item) => item.anchor.recordKey === recordKey)
      .map((item) => serializedCandidate(item, request));
    const choice = request.choices.find((item) => item.recordKey === recordKey);
    if (choice) {
      const selected = candidates.find(({ anchor: candidate }) => (
        JSON.stringify(candidate) === JSON.stringify(choice.anchor)
      ));
      return {
        kind: 'decision', recordKey,
        status: choice.kind === 'select' ? 'aligned' : 'skipped',
        rationale: choice.kind === 'select' ? 'user_selected' : 'skipped_by_user',
        reviewReason: null, anchor: selected?.anchor || null,
        orientationPolicy: choice.orientationPolicy,
        effectiveReverseComplement: selected
          ? selected.orientation[choice.orientationPolicy].effectiveReverseComplement : null,
        candidates
      };
    }
    if (!candidates.length || missing && recordKey === 'c') return {
      kind: 'decision', recordKey, status: 'skipped',
      rationale: 'skipped_no_candidate', reviewReason: null, anchor: null,
      orientationPolicy: 'preserve', effectiveReverseComplement: null, candidates: []
    };
    if (ambiguous && recordKey === 'b') return {
      kind: 'ambiguous', recordKey, candidates, directRbhCandidates: [],
      recommendedAnchor: candidates[0].anchor,
      recommendationReason: 'unique_representative'
    };
    return {
      kind: 'decision', recordKey, status: 'aligned',
      rationale: 'only_usable_candidate', reviewReason: 'only_usable_candidate',
      anchor: candidates[0].anchor, orientationPolicy: 'preserve',
      effectiveReverseComplement: false, candidates
    };
  });
  const status = records.some(({ kind }) => kind === 'ambiguous') ? 'ambiguous' : 'resolved';
  return {
    schema: 2, status, groupId: request.groupId, reference: request.reference,
    referenceDisplayedStrand: 1, referenceDisplayCenter: 65, records,
    plan: status === 'resolved' ? {
      schema: 2, groupId: request.groupId, reference: request.reference,
      records: records.map(({ recordKey, status: decisionStatus, rationale, anchor: chosen,
        orientationPolicy, effectiveReverseComplement }) => ({
        recordKey, status: decisionStatus, rationale, anchor: chosen,
        orientationPolicy, effectiveReverseComplement
      }))
    } : null
  };
};
const create = ({
  members = [reference, targetA, targetC], helper = (_operation, { request }) => ({
    result: responseFor(request)
  }), generation = async () => ({ status: 'ok' }), onError = null,
  previewCandidate = null, clearCandidatePreview = null,
  getCommittedRequest = () => renderRequest
} = {}) => {
  const controllerState = state();
  const currentGroup = { id: 'og-1', members, orthologEdges: [] };
  const helperCalls = [];
  const generationCalls = [];
  let svg = null;
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    getOrthogroupById: (id) => id === currentGroup.id ? currentGroup : null,
    getEnrichedOrthogroupMembers: () => currentGroup.members,
    getCommittedRequest,
    getCurrentSvg: () => svg,
    runAnalysis: async (options) => {
      generationCalls.push(structuredClone(options));
      const outcome = await generation(options);
      if (outcome.status === 'ok') {
        controllerState.similarityAlignmentPlan.value = options.canonicalStateOverride.similarityAlignmentPlan;
        controllerState.results.value = [{ name: 'aligned.svg' }];
        svg = null;
      }
      return outcome;
    },
    cancelRunAnalysis: () => {},
    runHelperOperation: async (operation, payload) => {
      helperCalls.push({ operation, payload: structuredClone(payload) });
      return helper(operation, payload, helperCalls.length);
    },
    resolveOperation: 'resolveSimilarityAlignment', onError,
    previewCandidate, clearCandidatePreview
  });
  return { actions, state: controllerState, helperCalls, generationCalls, currentGroup };
};
const startReview = (fixture, selected = reference) => fixture.actions.startFromPopup({
  groupId: 'og-1', reference: selected, mode: 'review'
});
const startAlign = (fixture, selected = reference) => fixture.actions.startFromPopup({
  groupId: 'og-1', reference: selected
});

test('normal Align applies only-usable and unusable targets with one helper and one Result', async () => {
  const applying = deferred();
  const fixture = create({ members: [reference, targetA, targetC],
    helper: (_operation, { request }) => {
      const response = responseFor(request);
      const unusable = response.records[2];
      unusable.status = 'skipped';
      unusable.rationale = 'skipped_unmappable';
      unusable.reviewReason = null;
      unusable.anchor = null;
      unusable.effectiveReverseComplement = null;
      unusable.candidates[0].usable = false;
      unusable.candidates[0].displayCenter = null;
      Object.assign(response.plan.records[2], {
        status: 'skipped', rationale: 'skipped_unmappable', anchor: null,
        effectiveReverseComplement: null
      });
      return { result: response };
    },
    generation: () => applying.promise });
  const priorResult = fixture.state.results.value;
  const operation = startAlign(fixture);
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(fixture.actions.status.value, 'applying');
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.actions.summary.value, null);
  assert.equal(fixture.state.results.value, priorResult);
  assert.equal(fixture.helperCalls.length, 1);
  assert.deepEqual(await startAlign(fixture), { status: 'busy' });
  applying.resolve({ status: 'ok' });
  assert.deepEqual(await operation, { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.generationCalls.length, 1);
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[1].rationale,
    'only_usable_candidate');
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[2].rationale,
    'skipped_unmappable');
  assert.deepEqual(fixture.actions.summary.value && [
    fixture.actions.summary.value.aligned, fixture.actions.summary.value.noCandidate,
    fixture.actions.summary.value.reversed
  ], [1, 1, 0]);
  assert.notEqual(fixture.state.results.value, priorResult);
});

test('normal Align keeps a missing target unchanged in the resolved plan', async () => {
  const fixture = create({ members: [reference, targetA],
    helper: (_operation, { request }) => ({ result: responseFor(request, { missing: true }) }) });
  assert.deepEqual(await startAlign(fixture), { status: 'ok' });
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[2].rationale,
    'skipped_no_candidate');
  assert.equal(fixture.actions.summary.value.noCandidate, 1);
});

test('normal Align uses Python unique-direct-RBH plan without a review or second helper', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => {
      assert.equal(request.members.filter(({ anchor: item }) => item.recordKey === 'b').length, 2);
      assert.equal(request.directEdges.length, 1);
      const response = responseFor(request);
      response.records[1].rationale = 'unique_direct_rbh';
      response.records[1].reviewReason = 'unique_direct_rbh';
      response.records[1].candidates[0].directEvidence = ['rbh'];
      response.plan.records[1].rationale = 'unique_direct_rbh';
      return { result: response };
    } });
  fixture.currentGroup.orthologEdges.push({
    orthogroupId: 'og-1', queryRecordIndex: 0, queryProteinId: 'clicked',
    subjectRecordIndex: 1, subjectProteinId: 'target-a', edgeKind: 'rbh'
  });
  assert.deepEqual(await startAlign(fixture), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.generationCalls.length, 1);
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[1].rationale,
    'unique_direct_rbh');
});

test('normal Align opens complete review for genuine Python ambiguity', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  assert.deepEqual(await startAlign(fixture), { status: 'reviewing' });
  assert.equal(fixture.generationCalls.length, 0);
  assert.equal(fixture.actions.draft.value.rows[0].reason, 'unique_representative');
  assert.equal(fixture.actions.draft.value.rows[0].choice.candidateKey,
    JSON.stringify(anchor('b', 'target-a', 1)));
  assert.equal(fixture.actions.canApply.value, true);
});

test('explicit review of a resolved response keeps the full local draft without generation', async () => {
  const fixture = create();
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
  assert.equal(fixture.actions.draft.value.response.status, 'resolved');
  assert.equal(fixture.actions.draft.value.rows[0].reason, 'only_usable_candidate');
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.generationCalls.length, 0);
});

test('automatic generation failure opens the resolved draft for one Apply retry', async () => {
  let fail = true;
  const fixture = create({ generation: async () => fail
    ? { status: 'error', error: new Error('Render failed.') }
    : { status: 'ok' } });
  const previous = fixture.state.results.value;
  assert.deepEqual(await startAlign(fixture), { status: 'error' });
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.dialogOpen.value, true);
  assert.equal(fixture.actions.draft.value.response.status, 'resolved');
  assert.match(fixture.actions.error.value.message, /Render failed/);
  assert.equal(fixture.state.results.value, previous);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  fail = false;
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 2);
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.actions.summary.value.aligned, 2);
});

test('canceled or stale resolved completion cannot install a plan or Result', async () => {
  const pendingHelper = deferred();
  const resolving = create({ helper: (_operation, { request }) => pendingHelper.promise.then(
    () => ({ result: responseFor(request) })
  ) });
  const first = startAlign(resolving);
  resolving.actions.cancel();
  pendingHelper.resolve();
  assert.deepEqual(await first, { status: 'stale' });
  assert.equal(resolving.generationCalls.length, 0);
  const pendingGeneration = deferred();
  const applying = create({ generation: () => pendingGeneration.promise });
  const prior = applying.state.results.value;
  const second = startAlign(applying);
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(applying.actions.status.value, 'applying');
  applying.actions.cancel();
  pendingGeneration.resolve({ status: 'canceled' });
  assert.deepEqual(await second, { status: 'stale' });
  assert.equal(applying.state.results.value, prior);
  assert.equal(applying.state.similarityAlignmentPlan.value, null);
  assert.equal(applying.actions.summary.value, null);
});

test('busy is synchronous, duplicate starts are ignored, and exact popup reference survives', async () => {
  const pending = deferred();
  const fixture = create({ members: [reference, otherReference, targetA, targetC],
    helper: (_operation, { request }) => pending.promise.then(() => ({ result: responseFor(request) })) });
  const operation = startReview(fixture);
  assert.equal(fixture.actions.status.value, 'resolving');
  assert.equal(fixture.actions.busy.value, true);
  assert.deepEqual(await startReview(fixture, otherReference), { status: 'busy' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.deepEqual(fixture.helperCalls[0].payload.request.reference,
    anchor('a', 'clicked', 4));
  pending.resolve();
  assert.deepEqual(await operation, { status: 'reviewing' });
  assert.equal(fixture.actions.busy.value, false);
  assert.equal(fixture.generationCalls.length, 0);
});

test('drawer preserves its selected exact member and opens the same review', async () => {
  const fixture = create({ members: [reference, otherReference, targetA, targetC] });
  const option = fixture.actions.drawerReferenceOptions('og-1').find(({ anchor: a }) => (
    a.biologicalFeatureId === 'other'
  ));
  fixture.actions.setDrawerReference('og-1', option.key);
  assert.deepEqual(await fixture.actions.startFromDrawer({ groupId: 'og-1', mode: 'review' }),
    { status: 'reviewing' });
  assert.deepEqual(fixture.helperCalls[0].payload.request.reference,
    anchor('a', 'other', 7));
});

test('explicit review opens an applicable Python-selected draft for either response status', async () => {
  for (const ambiguous of [false, true]) {
    const fixture = create({ members: [reference, targetA, targetB, targetC],
      helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous }) }) });
    assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
    assert.equal(fixture.actions.dialogOpen.value, true);
    assert.equal(fixture.actions.canApply.value, true);
    assert.equal(fixture.actions.unresolvedCount.value, 0);
    assert.deepEqual(fixture.actions.draft.value.rows.map(({ recordKey }) => recordKey), ['b', 'c']);
    assert.equal(fixture.actions.draft.value.rows[0].choice.candidateKey,
      JSON.stringify(anchor('b', 'target-a', 1)));
    assert.equal(fixture.actions.draft.value.rows[0].reason,
      ambiguous ? 'unique_representative' : 'only_usable_candidate');
    assert.equal(fixture.actions.draft.value.rows[0].orientationPolicy, 'preserve');
    assert.equal(fixture.actions.draft.value.reference.featureIdentifier, 'clicked');
    assert.equal(fixture.actions.draft.value.reference.coordinates, '51..70 bp');
    assert.equal(fixture.actions.draft.value.rows[0].recommendedKey,
      fixture.actions.draft.value.rows[0].choice.candidateKey);
    assert.match(fixture.actions.draft.value.rows[0].recommendationReasonLabel, /candidate/);
    assert.equal(fixture.generationCalls.length, 0);
  }
});

test('missing target is explicitly unchanged and complete Apply includes every target', async () => {
  const fixture = create({ members: [reference, targetA],
    helper: (_operation, { request }) => ({ result: responseFor(request, { missing: true }) }) });
  await startReview(fixture);
  assert.equal(fixture.actions.draft.value.rows[1].unchanged, true);
  assert.equal(fixture.actions.draft.value.rows[1].reason, 'skipped_no_candidate');
  assert.deepEqual(fixture.actions.skipRecord('c'), { status: 'rejected' });
  assert.equal(fixture.actions.draft.value.rows[1].reason, 'skipped_no_candidate');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 1);
  assert.deepEqual(fixture.helperCalls[1].payload.request.choices.map(({ recordKey, kind }) =>
    [recordKey, kind]), [['b', 'select'], ['c', 'skip']]);
  assert.equal(fixture.state.similarityAlignmentPlan.value.schema, 2);
  assert.deepEqual(fixture.state.results.value, [{ name: 'aligned.svg' }]);
});

test('candidate, canvas, Skip, and orientation use local edits with zero Worker jobs', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  await startReview(fixture);
  const count = fixture.helperCalls.length;
  assert.deepEqual(fixture.actions.selectCandidate('b', anchor('b', 'target-b', 2)),
    { status: 'selected' });
  assert.equal(fixture.actions.draft.value.rows[0].choice.candidateKey,
    JSON.stringify(anchor('b', 'target-b', 2)));
  assert.deepEqual(fixture.actions.setOrientation('b', 'match_reference'),
    { status: 'selected' });
  assert.equal(fixture.actions.draft.value.rows[0].orientationEffect, 'preserve');
  fixture.actions.skipRecord('b');
  assert.equal(fixture.actions.draft.value.rows[0].orientationPolicy, 'preserve');
  fixture.actions.selectCandidate('b', anchor('b', 'target-a', 1));
  fixture.actions.setOrientation('b', 'match_reference');
  assert.equal(fixture.actions.draft.value.rows[0].orientationEffect, 'reverse_whole_record');
  fixture.actions.setOrientation('c', 'match_reference');
  assert.equal(fixture.actions.draft.value.rows[1].orientationEffect, 'preserve_unknown_strand');
  assert.equal(fixture.helperCalls.length, count);
  assert.equal(fixture.generationCalls.length, 0);
});

test('Apply validates one explicit batch and commits through one generation call', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  await startReview(fixture);
  fixture.actions.setOrientation('b', 'match_reference');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 1);
  assert.deepEqual(fixture.helperCalls[1].payload.request.choices, [
    { recordKey: 'b', kind: 'select', anchor: anchor('b', 'target-a', 1),
      orientationPolicy: 'match_reference' },
    { recordKey: 'c', kind: 'select', anchor: anchor('c', 'target-c', 3),
      orientationPolicy: 'preserve' }
  ]);
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[1].effectiveReverseComplement, true);
  assert.equal(fixture.actions.status.value, 'idle');
});

test('duplicate Apply while validation is pending starts one batch', async () => {
  const pending = deferred();
  const fixture = create({ helper: (_operation, { request }, count) => (
    count === 2 ? pending.promise.then(() => ({ result: responseFor(request) }))
      : { result: responseFor(request) }
  ) });
  await startReview(fixture);
  const operation = fixture.actions.applyDraft();
  assert.equal(fixture.actions.status.value, 'applying');
  assert.equal(fixture.actions.busy.value, true);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'rejected' });
  assert.deepEqual(await startReview(fixture), { status: 'busy' });
  pending.resolve();
  assert.deepEqual(await operation, { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 1);
});

test('validation and generation failures preserve editable draft and prior artifact', async () => {
  for (const failure of ['validation', 'generation']) {
    let fail = true;
    const fixture = create({
      helper: (_operation, { request }, count) => {
        if (failure === 'validation' && count > 1 && fail) throw new Error('Invalid choice. Select another.');
        return { result: responseFor(request) };
      },
      generation: async () => failure === 'generation' && fail
        ? { status: 'error', error: new Error('Generation failed. Retry Apply.') }
        : { status: 'ok' }
    });
    await startReview(fixture);
    const before = fixture.actions.draft.value;
    assert.equal((await fixture.actions.applyDraft()).status, 'error');
    assert.equal(fixture.actions.status.value, 'reviewing');
    assert.equal(fixture.actions.draft.value, before);
    assert.deepEqual(fixture.state.results.value, [{ name: 'prior.svg' }]);
    assert.equal(fixture.state.similarityAlignmentPlan.value, null);
    assert.match(fixture.actions.error.value.message, /Select another|Retry Apply/);
    fail = false;
    fixture.actions.skipRecord('b');
    assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  }
});

test('Cancel preserves committed diagram, plan, and Result', async () => {
  const fixture = create();
  await startReview(fixture);
  const previous = fixture.state.results.value;
  assert.deepEqual(fixture.actions.cancel(), { status: 'canceled' });
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.equal(fixture.state.results.value, previous);
  assert.equal(fixture.generationCalls.length, 0);
});

test('stale initial resolution and canceled Apply cannot publish late completions', async () => {
  const initial = deferred();
  const resolving = create({ helper: (_operation, { request }) => initial.promise.then(() => ({
    result: responseFor(request)
  })) });
  const first = startReview(resolving);
  resolving.actions.cancel();
  initial.resolve();
  assert.deepEqual(await first, { status: 'stale' });
  assert.equal(resolving.actions.draft.value, null);
  const apply = deferred();
  const fixture = create({ helper: (_operation, { request }, count) => count === 2
    ? apply.promise.then(() => ({ result: responseFor(request) }))
    : { result: responseFor(request) } });
  await startReview(fixture);
  const second = fixture.actions.applyDraft();
  fixture.actions.cancel();
  apply.resolve();
  assert.deepEqual(await second, { status: 'stale' });
  assert.equal(fixture.generationCalls.length, 0);
  assert.deepEqual(fixture.state.results.value, [{ name: 'prior.svg' }]);
});

test('changed committed Result rejects an old draft before validation', async () => {
  const fixture = create();
  await startReview(fixture);
  fixture.state.results.value = [{ name: 'new.svg' }];
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'stale' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.deepEqual(fixture.state.results.value, [{ name: 'new.svg' }]);
});

test('active schema-2 plan validates on regeneration and supports Reset and stable reorder', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.activePlanInspector.value.modeLabel, 'Align');
  const calls = fixture.helperCalls.length;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, calls + 1);
  fixture.actions.retainForStableReorder(['c', 'a', 'b']);
  assert.deepEqual(fixture.state.similarityAlignmentPlan.value.records.map(({ recordKey }) => recordKey),
    ['c', 'a', 'b']);
  assert.deepEqual((await fixture.actions.resetAlignment()).status, 'ok');
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
});

test('Python recommendation and direct-RBH reason control initial selection', async () => {
  const ambiguous = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => {
      const response = responseFor(request, { ambiguous: true });
      response.records[1].recommendedAnchor = anchor('b', 'target-b', 2);
      response.records[1].recommendationReason = 'deterministic_candidate_1';
      return { result: response };
    } });
  await startReview(ambiguous);
  assert.equal(ambiguous.actions.draft.value.rows[0].choice.candidateKey,
    JSON.stringify(anchor('b', 'target-b', 2)));
  assert.equal(ambiguous.actions.draft.value.rows[0].reason, 'deterministic_candidate_1');
  const direct = create({ helper: (_operation, { request }) => {
    const response = responseFor(request);
    response.records[1].rationale = 'unique_direct_rbh';
    response.records[1].reviewReason = 'unique_direct_rbh';
    response.plan.records[1].rationale = 'unique_direct_rbh';
    return { result: response };
  } });
  await startReview(direct);
  assert.equal(direct.actions.draft.value.rows[0].reason, 'unique_direct_rbh');
});

test('stale active target enters the same repair draft and needs an explicit new choice', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const priorResult = fixture.state.results.value;
  fixture.currentGroup.members.splice(1, 1, targetB);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-target' });
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.unresolvedCount.value, 1);
  assert.equal(fixture.actions.canApply.value, false);
  assert.equal(fixture.state.results.value, priorResult);
  assert.equal(fixture.actions.selectCandidate('b', anchor('b', 'target-b', 2)).status,
    'selected');
  assert.equal(fixture.actions.canApply.value, true);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[1].anchor.biologicalFeatureId,
    'target-b');
});

test('drawer requires a unique exact reference before starting the Worker', async () => {
  const fixture = create();
  assert.match(fixture.actions.drawerDisabledReason('og-1'), /Select an exact reference/);
  assert.deepEqual(await fixture.actions.startFromDrawer({ groupId: 'og-1' }),
    { status: 'rejected' });
  assert.equal(fixture.helperCalls.length, 0);
  const duplicate = create({ members: [reference, { ...reference }, targetA] });
  assert.deepEqual(duplicate.actions.drawerReferenceOptions('og-1')
    .filter(({ anchor: selected }) => selected.biologicalFeatureId === 'clicked'), []);
});

test('malformed Python projection and initial Worker errors leave the prior Result intact', async () => {
  for (const helper of [
    (_operation, { request }) => ({ result: { ...responseFor(request), schema: 1 } }),
    () => { throw new Error('Python resolver unavailable. Retry Align.'); }
  ]) {
    const fixture = create({ helper });
    const prior = fixture.state.results.value;
    assert.equal((await startReview(fixture)).status, 'error');
    assert.equal(fixture.actions.status.value, 'idle');
    assert.equal(fixture.actions.draft.value, null);
    assert.equal(fixture.state.results.value, prior);
    assert.equal(fixture.generationCalls.length, 0);
    assert.ok(fixture.actions.error.value);
  }
});

test('candidate preview delegates to the existing owner and Cancel clears it', async () => {
  const previewed = [];
  let cleared = 0;
  const fixture = create({ previewCandidate: (value) => previewed.push(value),
    clearCandidatePreview: () => { cleared += 1; } });
  await startReview(fixture);
  fixture.actions.previewCandidate(anchor('b', 'target-a', 1));
  assert.deepEqual(previewed, [anchor('b', 'target-a', 1)]);
  fixture.actions.cancel();
  assert.ok(cleared >= 1);
});

test('manual orientation and invalidating edits clear an active plan with a notice', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.setManualOrientation(fixture.state.linearSeqs[1], true), true);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.equal(fixture.state.linearSeqs[1].region_reverse, true);
  assert.match(fixture.actions.notice.value, /orientation changed/);
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.clearForMutation('record crop changed.'), true);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.match(fixture.actions.notice.value, /record crop changed/);
});

test('missing saved reference blocks regeneration while preserving plan and Result', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const previous = fixture.state.results.value;
  const plan = fixture.state.similarityAlignmentPlan.value;
  fixture.currentGroup.members.splice(0, 1);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.state.results.value, previous);
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  assert.equal(fixture.actions.repair.value.kind, 'reference');
});

test('a committed schema-2 plan regenerates after its transient groups leave the catalog', async () => {
  let committed = renderRequest;
  const fixture = create({ getCommittedRequest: () => committed });
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const plan = fixture.state.similarityAlignmentPlan.value;
  committed = { ...renderRequest, layout: { similarityAlignment: structuredClone(plan) } };
  fixture.currentGroup.id = 'no-longer-in-catalog';
  const calls = fixture.helperCalls.length;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, calls);
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  fixture.actions.retainForStableReorder(['c', 'a', 'b']);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  committed.layout.similarityAlignment.records[1].anchor.biologicalFeatureId = 'different';
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
});

test('candidate labels retain biological metadata and displayed record names', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  fixture.state.linearSeqs[1].definition = 'Genome B';
  await startReview(fixture);
  const row = fixture.actions.draft.value.rows[0];
  assert.equal(row.recordLabel, 'Genome B');
  assert.equal(row.candidates[0].label, 'gene-target-a');
  assert.equal(row.candidates[0].coordinates, '21..40 bp');
  assert.equal(row.candidates[0].displayedStrand, '-');
});
