import assert from 'node:assert/strict';
import test from 'node:test';

const ref = (value) => ({ value });
globalThis.window = {
  Vue: {
    ref,
    computed: (getter) => ({ get value() { return getter(); } })
  }
};

const { createSimilarityAlignmentActions } = await import(
  '../../gbdraw/web/js/app/similarity-alignment.js'
);

const anchor = (recordKey, biologicalFeatureId, sourceFeatureIndex) => ({
  recordKey,
  biologicalFeatureId,
  sourceFeatureIndex,
  stableFeatureSvgId: `stable-${biologicalFeatureId}`
});

const reference = {
  ...anchor('record-a', 'clicked-inparalog', 4),
  recordIndex: 0,
  featureIndex: 4,
  stableFeatureSvgId: 'stable-clicked-inparalog',
  start: 10,
  end: 40,
  strand: '+',
  proteinId: 'protein-a',
  sourceProteinId: 'protein-a',
  representative: false,
  role: 'inparalog'
};
const targetA = {
  ...anchor('record-b', 'target-a', 1),
  recordIndex: 1,
  featureIndex: 1,
  stableFeatureSvgId: 'stable-target-a',
  start: 100,
  end: 130,
  strand: '-',
  proteinId: 'protein-b1',
  sourceProteinId: 'protein-b1',
  representative: true,
  role: 'anchor'
};
const targetB = {
  ...anchor('record-b', 'target-b', 2),
  recordIndex: 1,
  featureIndex: 2,
  stableFeatureSvgId: 'stable-target-b',
  start: 160,
  end: 190,
  strand: '+',
  proteinId: 'protein-b2',
  sourceProteinId: 'protein-b2',
  representative: false,
  role: 'inparalog'
};

const renderRequest = {
  mode: 'linear',
  records: [
    { recordKey: 'record-a', region: null, presentation: { reverseComplement: false } },
    { recordKey: 'record-b', region: null, presentation: { reverseComplement: false } }
  ]
};

const state = () => ({
  featureCatalog: ref({
    items: [{
      recordKeys: ['record-a', 'record-b'],
      sequenceSources: [{ sequence: 'A'.repeat(300) }, { sequence: 'C'.repeat(300) }]
    }]
  }),
  similarityAlignmentPlan: ref(null),
  linearRecordTranslations: ref([]),
  results: ref([{ name: 'committed.svg' }])
});

const decision = (recordKey, status, rationale, selectedAnchor = null) => ({
  kind: 'decision',
  recordKey,
  status,
  rationale,
  anchor: selectedAnchor,
  effectiveReverseComplement: null
});

const planDecision = ({ kind: _kind, ...value }) => value;

const resolvedResponse = (request, selected = anchor('record-b', 'target-a', 1)) => {
  const records = [
    decision('record-a', 'reference', 'reference', request.reference),
    decision('record-b', 'aligned', 'only_usable_candidate', selected)
  ];
  return {
    schema: 1,
    status: 'resolved',
    mode: request.mode,
    groupId: request.groupId,
    reference: request.reference,
    records,
    plan: {
      schema: 1,
      mode: request.mode,
      groupId: request.groupId,
      reference: request.reference,
      records: records.map(planDecision)
    }
  };
};

const ambiguityResponse = (request) => ({
  schema: 1,
  status: 'ambiguous',
  mode: request.mode,
  groupId: request.groupId,
  reference: request.reference,
  records: [
    decision('record-a', 'reference', 'reference', request.reference),
    {
      kind: 'ambiguous',
      recordKey: 'record-b',
      candidates: [targetA, targetB].map((member) => ({
        anchor: anchor(member.recordKey, member.biologicalFeatureId, member.sourceFeatureIndex),
        displayedStrand: member.strand === '-' ? -1 : 1,
        hidden: false,
        representative: member.representative,
        role: member.role
      })),
      directRbhCandidates: []
    }
  ],
  plan: null
});

const group = (members = [reference, targetA], orthologEdges = [], relatedEdges = []) => ({
  id: 'og-1',
  members,
  orthologEdges,
  relatedEdges
});

const create = ({
  members = [reference, targetA],
  orthologEdges = [],
  relatedEdges = [],
  helper,
  runAnalysis = async () => ({ status: 'ok' }),
  cancelRunAnalysis = () => {}
} = {}) => {
  const controllerState = state();
  const currentGroup = group(members, orthologEdges, relatedEdges);
  const helperCalls = [];
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    getOrthogroupById: (id) => id === currentGroup.id ? currentGroup : null,
    getEnrichedOrthogroupMembers: () => currentGroup.members,
    getCommittedRequest: () => renderRequest,
    runAnalysis,
    cancelRunAnalysis,
    resolveOperation: 'resolveSimilarityAlignment',
    runHelperOperation: async (operation, payload) => {
      helperCalls.push({ operation, payload: structuredClone(payload) });
      return helper(operation, payload, helperCalls.length);
    }
  });
  return { actions, state: controllerState, helperCalls };
};

test('popup sends its exact clicked reference and auto-resolved Apply uses one generation transaction', async () => {
  let generationCalls = 0;
  let observedOverride = null;
  const fixture = create({
    orthologEdges: [{
      orthogroupId: 'og-1',
      queryProteinId: 'protein-a',
      subjectProteinId: 'protein-b1',
      queryRecordIndex: 0,
      subjectRecordIndex: 1,
      edgeKind: 'rbh'
    }],
    relatedEdges: [{ edgeKind: 'path-only-evidence-must-not-cross-helper' }],
    helper: async (_operation, { request }) => ({ result: resolvedResponse(request) }),
    runAnalysis: async ({ canonicalStateOverride }) => {
      generationCalls += 1;
      observedOverride = canonicalStateOverride;
      assert.equal(fixture.state.similarityAlignmentPlan.value, null);
      assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
      return { status: 'ok' };
    }
  });

  assert.deepEqual(await fixture.actions.startFromPopup({
    groupId: 'og-1', reference, mode: 'position'
  }), { status: 'ok' });
  assert.equal(generationCalls, 1);
  assert.equal(fixture.helperCalls[0].operation, 'resolveSimilarityAlignment');
  assert.deepEqual(fixture.helperCalls[0].payload.request.reference, anchor(
    'record-a', 'clicked-inparalog', 4
  ));
  assert.deepEqual(fixture.helperCalls[0].payload.request.directEdges, [{
    groupId: 'og-1',
    query: anchor('record-a', 'clicked-inparalog', 4),
    subject: anchor('record-b', 'target-a', 1),
    edgeKind: 'rbh'
  }]);
  assert.deepEqual(observedOverride.linearRecordTranslations, [
    { recordKey: 'record-a', x: 0, y: 0 },
    { recordKey: 'record-b', x: 0, y: 0 }
  ]);
  assert.equal(observedOverride.similarityAlignmentPlan.groupId, 'og-1');
});

test('drawer rejects a group-only action before helper or generation', async () => {
  const fixture = create({
    helper: async () => { throw new Error('must not run'); },
    runAnalysis: async () => { throw new Error('must not generate'); }
  });

  assert.deepEqual(await fixture.actions.startFromDrawer({ groupId: 'og-1' }), {
    status: 'rejected'
  });
  assert.equal(fixture.helperCalls.length, 0);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
});

test('ambiguity is immutable and Select returns to the same resolver before explicit Apply', async () => {
  let generationCalls = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }, call) => {
      if (call === 1) return { result: ambiguityResponse(request) };
      const selected = request.choices[0].anchor;
      const response = resolvedResponse(request, selected);
      response.records[1].rationale = 'user_selected';
      response.plan.records[1].rationale = 'user_selected';
      return { result: response };
    },
    runAnalysis: async () => {
      generationCalls += 1;
      return { status: 'ok' };
    }
  });

  assert.deepEqual(await fixture.actions.startFromPopup({
    groupId: 'og-1', reference, mode: 'position'
  }), { status: 'ambiguous' });
  assert.equal(Object.isFrozen(fixture.actions.draft.value), true);
  assert.equal(generationCalls, 0);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);

  const selected = anchor('record-b', 'target-b', 2);
  assert.deepEqual(await fixture.actions.selectCandidate('record-b', selected), {
    status: 'ready'
  });
  assert.equal(fixture.helperCalls.length, 2);
  assert.deepEqual(fixture.helperCalls[1].payload.request.choices, [{
    recordKey: 'record-b', kind: 'select', anchor: selected
  }]);
  assert.equal(generationCalls, 0);
  assert.equal(fixture.actions.canApply.value, true);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(generationCalls, 1);
});

test('Skip returns to the Python resolver and does not synthesize a plan in JavaScript', async () => {
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }, call) => {
      if (call === 1) return { result: ambiguityResponse(request) };
      assert.deepEqual(request.choices, [{ recordKey: 'record-b', kind: 'skip', anchor: null }]);
      const records = [
        decision('record-a', 'reference', 'reference', request.reference),
        decision('record-b', 'skipped', 'skipped_by_user', null)
      ];
      return { result: {
        schema: 1, status: 'resolved', mode: request.mode, groupId: request.groupId,
        reference: request.reference, records,
        plan: { schema: 1, mode: request.mode, groupId: request.groupId,
          reference: request.reference, records: records.map(planDecision) }
      } };
    }
  });
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  assert.deepEqual(await fixture.actions.skipRecord('record-b'), { status: 'ready' });
  assert.equal(fixture.actions.draft.value.response.plan.records[1].rationale, 'skipped_by_user');
});

test('Cancel and stale helper completion discard drafts without changing the committed artifact', async () => {
  let resolveHelper;
  let generationCalls = 0;
  const fixture = create({
    helper: async (_operation, { request }) => new Promise((resolve) => {
      resolveHelper = () => resolve({ result: ambiguityResponse(request) });
    }),
    runAnalysis: async () => { generationCalls += 1; return { status: 'ok' }; }
  });
  const pending = fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  await Promise.resolve();
  assert.deepEqual(fixture.actions.cancel(), { status: 'canceled' });
  resolveHelper();
  assert.deepEqual(await pending, { status: 'stale' });
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
  assert.equal(generationCalls, 0);
});

for (const failedStatus of ['error', 'canceled']) {
  test(`render ${failedStatus} discards the draft and preserves canonical state`, async () => {
    const fixture = create({
      helper: async (_operation, { request }) => ({ result: resolvedResponse(request) }),
      runAnalysis: async () => ({ status: failedStatus })
    });
    assert.deepEqual(await fixture.actions.startFromPopup({ groupId: 'og-1', reference }), {
      status: failedStatus
    });
    assert.equal(fixture.state.similarityAlignmentPlan.value, null);
    assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
    assert.equal(fixture.actions.draft.value, null);
  });
}

test('unknown helper response values are rejected without generation', async () => {
  let generationCalls = 0;
  const fixture = create({
    helper: async (_operation, { request }) => {
      const response = resolvedResponse(request);
      response.records[1].status = 'ranked';
      return { result: response };
    },
    runAnalysis: async () => { generationCalls += 1; return { status: 'ok' }; }
  });
  assert.deepEqual(await fixture.actions.startFromPopup({ groupId: 'og-1', reference }), {
    status: 'error'
  });
  assert.equal(generationCalls, 0);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
});

test('Worker error discards the draft and preserves the committed artifact', async () => {
  let generationCalls = 0;
  const fixture = create({
    helper: async () => { throw new Error('injected alignment helper failure'); },
    runAnalysis: async () => { generationCalls += 1; return { status: 'ok' }; }
  });
  assert.deepEqual(await fixture.actions.startFromPopup({ groupId: 'og-1', reference }), {
    status: 'error'
  });
  assert.match(fixture.actions.error.value.message, /injected alignment helper failure/);
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(generationCalls, 0);
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
});

test('a newer action cancels an in-flight Apply before publishing its own draft', async () => {
  let finishApply;
  let cancelCalls = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }, call) => ({
      result: call === 1 ? resolvedResponse(request) : ambiguityResponse(request)
    }),
    runAnalysis: async () => new Promise((resolve) => { finishApply = resolve; }),
    cancelRunAnalysis: () => { cancelCalls += 1; }
  });

  const first = fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  while (typeof finishApply !== 'function') await Promise.resolve();
  const second = fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  assert.equal(cancelCalls, 1);
  finishApply({ status: 'canceled' });
  assert.deepEqual(await first, { status: 'stale' });
  assert.deepEqual(await second, { status: 'ambiguous' });
  assert.equal(fixture.actions.status.value, 'ambiguous');
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
});
