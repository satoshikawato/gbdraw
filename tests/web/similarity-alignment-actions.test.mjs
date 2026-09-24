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
const targetC = {
  ...anchor('record-c', 'target-c', 3),
  recordIndex: 2,
  featureIndex: 3,
  stableFeatureSvgId: 'stable-target-c',
  start: 210,
  end: 250,
  strand: '-',
  proteinId: 'protein-c1',
  sourceProteinId: 'protein-c1',
  representative: false,
  role: 'coortholog'
};
const targetD = {
  ...anchor('record-c', 'target-d', 5),
  recordIndex: 2,
  featureIndex: 5,
  stableFeatureSvgId: 'stable-target-d',
  start: 260,
  end: 290,
  strand: '+',
  proteinId: 'protein-c2',
  sourceProteinId: 'protein-c2',
  representative: true,
  role: 'anchor'
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
  linearSeqs: [
    { uid: 'record-a', region_reverse: false },
    { uid: 'record-b', region_reverse: false }
  ],
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
  cancelRunAnalysis = () => {},
  previewCandidate = null,
  clearCandidatePreview = null,
  getCurrentSvg = null
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
    getCurrentSvg,
    previewCandidate,
    clearCandidatePreview,
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
  assert.equal(fixture.actions.dialogOpen.value, false);
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

test('drawer requires and resolves one exact reference before either explicit action', async () => {
  const fixture = create({
    helper: async (_operation, { request }) => ({ result: resolvedResponse(request) })
  });
  const options = fixture.actions.drawerReferenceOptions('og-1');
  assert.equal(options.length, 2);
  assert.match(options[0].label, /record-a.*clicked-inparalog.*11\.\.40.*\(\+\)/);
  assert.match(fixture.actions.drawerDisabledReason('og-1'), /Select an exact reference/);
  assert.equal(fixture.actions.setDrawerReference('og-1', options[0].key), true);
  assert.equal(fixture.actions.drawerDisabledReason('og-1'), '');
  assert.deepEqual(await fixture.actions.startFromDrawer({
    groupId: 'og-1', mode: 'position_and_orientation'
  }), { status: 'ok' });
  assert.equal(fixture.helperCalls[0].payload.request.mode, 'position_and_orientation');
  assert.deepEqual(fixture.helperCalls[0].payload.request.reference, options[0].anchor);
});

test('drawer does not offer duplicate canonical feature identities as exact references', () => {
  const fixture = create({
    members: [reference, {
      ...reference,
      sourceFeatureIndex: 5,
      featureIndex: 5,
      stableFeatureSvgId: 'stable-duplicate-clicked-inparalog'
    }, targetA],
    helper: async () => { throw new Error('must not resolve'); }
  });

  assert.deepEqual(
    fixture.actions.drawerReferenceOptions('og-1').map(({ anchor: value }) => value),
    [anchor('record-b', 'target-a', 1)]
  );
});

test('Select stays local until one explicit batch resolver Apply', async () => {
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
  assert.deepEqual(fixture.actions.selectCandidate('record-b', selected), { status: 'selected' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(generationCalls, 0);
  assert.equal(fixture.actions.canApply.value, true);
  assert.equal(fixture.actions.draft.value.response.plan, null);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.deepEqual(fixture.helperCalls[1].payload.request.choices, [{
    recordKey: 'record-b', kind: 'select', anchor: selected
  }]);
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(generationCalls, 1);
});

test('Skip remains local and Apply uses the Python plan', async () => {
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
  assert.deepEqual(fixture.actions.skipRecord('record-b'), { status: 'selected' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.actions.draft.value.response.plan, null);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
});

test('dialog state retains every resolver-reported ambiguity, facts, choices, and Apply reason', async () => {
  const threeRecordRequest = {
    ...renderRequest,
    records: [
      ...renderRequest.records,
      { recordKey: 'record-c', region: null, presentation: { reverseComplement: false } }
    ]
  };
  const currentGroup = group(
    [reference, targetA, targetB, targetC, targetD],
    [{
      orthogroupId: 'og-1', queryProteinId: 'protein-a', subjectProteinId: 'protein-b1',
      queryRecordIndex: 0, subjectRecordIndex: 1, edgeKind: 'coortholog'
    }]
  );
  const helperCalls = [];
  const controllerState = state();
  const responseFor = (request) => {
    const decisions = [decision('record-a', 'reference', 'reference', request.reference)];
    const records = [decisions[0]];
    for (const [recordKey, candidates] of [
      ['record-b', [targetA, targetB]], ['record-c', [targetC, targetD]]
    ]) {
      const choice = request.choices.find((entry) => entry.recordKey === recordKey);
      if (choice) {
        const resolved = choice.kind === 'skip'
          ? decision(recordKey, 'skipped', 'skipped_by_user', null)
          : decision(recordKey, 'aligned', 'user_selected', choice.anchor);
        decisions.push(resolved);
        records.push(resolved);
      } else {
        records.push({
          kind: 'ambiguous', recordKey,
          candidates: candidates.map((member) => ({
            anchor: anchor(member.recordKey, member.biologicalFeatureId, member.sourceFeatureIndex),
            displayedStrand: member.strand === '-' ? -1 : 1,
            hidden: false, representative: member.representative, role: member.role
          })),
          directRbhCandidates: []
        });
      }
    }
    return {
      schema: 1,
      status: records.some(({ kind }) => kind === 'ambiguous') ? 'ambiguous' : 'resolved',
      mode: request.mode,
      groupId: request.groupId,
      reference: request.reference,
      records,
      plan: records.some(({ kind }) => kind === 'ambiguous') ? null : {
        schema: 1, mode: request.mode, groupId: request.groupId,
        reference: request.reference, records: decisions.map(planDecision)
      }
    };
  };
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    getOrthogroupById: () => currentGroup,
    getEnrichedOrthogroupMembers: () => currentGroup.members,
    getCommittedRequest: () => threeRecordRequest,
    runAnalysis: async () => ({ status: 'ok' }),
    resolveOperation: 'resolveSimilarityAlignment',
    runHelperOperation: async (_operation, { request }) => {
      helperCalls.push(structuredClone(request));
      return { result: responseFor(request) };
    }
  });

  await actions.startFromPopup({ groupId: 'og-1', reference });
  assert.equal(actions.dialogOpen.value, true);
  assert.equal(actions.canApply.value, false);
  assert.match(actions.applyDisabledReason.value, /2 ambiguous records/);
  assert.deepEqual(actions.draft.value.ambiguities.map(({ recordKey }) => recordKey), [
    'record-b', 'record-c'
  ]);
  assert.deepEqual(actions.draft.value.ambiguities[0].candidates[0], {
    key: JSON.stringify(anchor('record-b', 'target-a', 1)),
    anchor: anchor('record-b', 'target-a', 1),
    featureIdentifier: 'target-a',
    label: 'CDS',
    coordinates: '101..130 bp',
    displayedStrand: '-',
    representative: true,
    role: 'anchor',
    directEvidence: ['COORTHOLOG']
  });

  assert.deepEqual(await actions.applyDraft(), { status: 'rejected' });
  assert.equal(helperCalls.length, 1);
  await actions.selectCandidate('record-b', anchor('record-b', 'target-b', 2));
  assert.equal(actions.draft.value.ambiguities.length, 2);
  assert.deepEqual(actions.draft.value.ambiguities[0].choice, {
    kind: 'select', candidateKey: JSON.stringify(anchor('record-b', 'target-b', 2))
  });
  assert.match(actions.applyDisabledReason.value, /1 ambiguous record/);
  await actions.selectCandidate('record-b', anchor('record-b', 'target-a', 1));
  await actions.selectCandidate('record-b', anchor('record-b', 'target-b', 2));
  await actions.skipRecord('record-c');
  await actions.selectCandidate('record-c', anchor('record-c', 'target-d', 5));
  await actions.skipRecord('record-c');
  assert.equal(helperCalls.length, 1);
  assert.equal(actions.canApply.value, true);
  assert.equal(actions.applyDisabledReason.value, '');
  assert.deepEqual(await actions.applyDraft(), { status: 'ok' });
  assert.equal(helperCalls.length, 2);
  assert.deepEqual(helperCalls.at(-1).choices.map(({ recordKey, kind }) => [recordKey, kind]), [
    ['record-b', 'select'], ['record-c', 'skip']
  ]);
});

test('successful Apply publishes the five-count live summary and source-relative plan inspector', async () => {
  let appliedPlan = null;
  const controllerState = state();
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    getOrthogroupById: () => group(),
    getEnrichedOrthogroupMembers: () => [reference, targetA],
    getCommittedRequest: () => ({
      ...renderRequest,
      records: [renderRequest.records[0], {
        ...renderRequest.records[1], presentation: { reverseComplement: false }
      }]
    }),
    runAnalysis: async ({ canonicalStateOverride }) => {
      appliedPlan = structuredClone(canonicalStateOverride.similarityAlignmentPlan);
      controllerState.similarityAlignmentPlan.value = structuredClone(appliedPlan);
      return { status: 'ok' };
    },
    resolveOperation: 'resolveSimilarityAlignment',
    runHelperOperation: async (_operation, { request }) => {
      const response = resolvedResponse(request);
      response.records[1].effectiveReverseComplement = true;
      response.plan.records[1].effectiveReverseComplement = true;
      return { result: response };
    }
  });
  assert.deepEqual(await actions.startFromPopup({
    groupId: 'og-1', reference, mode: 'position_and_orientation'
  }), { status: 'ok' });
  assert.deepEqual(actions.summary.value, {
    aligned: 1, unchanged: 1, explicitlySkipped: 0, noCandidate: 0, reversed: 1,
    text: 'Alignment applied: 1 aligned, 1 unchanged, 0 explicitly skipped, 0 no-candidate, 1 reversed.'
  });
  assert.equal(appliedPlan.mode, 'position_and_orientation');
  assert.equal(actions.activePlanInspector.value.reference.label, 'clicked-inparalog');
  assert.deepEqual(
    actions.activePlanInspector.value.records.map(({ anchorLabel, rationaleLabel, reversedFromSource }) => (
      [anchorLabel, rationaleLabel, reversedFromSource]
    )),
    [
      ['clicked-inparalog', 'Exact reference', false],
      ['target-a', 'Only usable candidate', true]
    ]
  );
});

test('plan inspector shows Skip rationale and derives rev from effective source orientation', () => {
  const controllerState = state();
  controllerState.similarityAlignmentPlan.value = {
    schema: 1,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4),
    records: [
      planDecision(decision(
        'record-a', 'reference', 'reference', anchor('record-a', 'clicked-inparalog', 4)
      )),
      planDecision(decision('record-b', 'skipped', 'skipped_by_user', null))
    ]
  };
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    getOrthogroupById: () => group(),
    getEnrichedOrthogroupMembers: () => [reference, targetA],
    getCommittedRequest: () => ({
      ...renderRequest,
      records: [renderRequest.records[0], {
        ...renderRequest.records[1], presentation: { reverseComplement: true }
      }]
    }),
    runAnalysis: async () => ({ status: 'ok' }),
    resolveOperation: 'resolveSimilarityAlignment',
    runHelperOperation: async () => { throw new Error('must not resolve'); }
  });

  assert.deepEqual(
    actions.activePlanInspector.value.records.map((record) => ({
      anchor: record.anchorLabel,
      rationale: record.rationaleLabel,
      rev: record.reversedFromSource
    })),
    [
      { anchor: 'clicked-inparalog', rationale: 'Exact reference', rev: false },
      { anchor: 'Skip', rationale: 'Skipped by user', rev: true }
    ]
  );
});

test('candidate preview delegates to the highlighting owner and Cancel restores it', async () => {
  const previewed = [];
  let clears = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }) => ({ result: ambiguityResponse(request) }),
    previewCandidate: (value) => previewed.push(value),
    clearCandidatePreview: () => { clears += 1; }
  });
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  fixture.actions.previewCandidate(anchor('record-b', 'target-b', 2));
  assert.deepEqual(previewed, [anchor('record-b', 'target-b', 2)]);
  fixture.actions.cancel();
  assert.ok(clears >= 2);
  assert.equal(fixture.actions.dialogOpen.value, false);
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

test('Reset removes only the active overlay and restores its immediate base values', async () => {
  let observedOverride = null;
  const fixture = create({
    helper: async () => { throw new Error('must not resolve'); },
    runAnalysis: async ({ canonicalStateOverride }) => {
      observedOverride = structuredClone(canonicalStateOverride);
      fixture.state.similarityAlignmentPlan.value = canonicalStateOverride.similarityAlignmentPlan;
      fixture.state.linearRecordTranslations.value = structuredClone(
        canonicalStateOverride.linearRecordTranslations
      );
      return { status: 'ok' };
    }
  });
  fixture.state.linearRecordTranslations.value = [
    { recordKey: 'record-a', x: 11, y: 2 },
    { recordKey: 'record-b', x: -7, y: 3 }
  ];
  fixture.state.similarityAlignmentPlan.value = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;

  assert.deepEqual(await fixture.actions.resetAlignment(), { status: 'ok' });
  assert.equal(observedOverride.similarityAlignmentPlan, null);
  assert.deepEqual(observedOverride.linearRecordTranslations, [
    { recordKey: 'record-a', x: 11, y: 2 },
    { recordKey: 'record-b', x: -7, y: 3 }
  ]);
  assert.match(fixture.actions.notice.value, /immediate pre-align baseline/);
  assert.deepEqual(await fixture.actions.resetAlignment(), { status: 'noop' });
});

test('manual orientation materializes effective orientation before clearing the plan', () => {
  const fixture = create({ helper: async () => { throw new Error('must not resolve'); } });
  const response = resolvedResponse({
    ...renderRequest,
    mode: 'position_and_orientation',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  });
  response.plan.records[1].effectiveReverseComplement = true;
  fixture.state.similarityAlignmentPlan.value = response.plan;

  assert.equal(fixture.actions.setManualOrientation(fixture.state.linearSeqs[1], false), true);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.equal(fixture.state.linearSeqs[1].region_reverse, false);
  assert.match(fixture.actions.notice.value, /record orientation changed/);
});

test('stable reorder remaps the plan and translations by record key', () => {
  const fixture = create({ helper: async () => { throw new Error('must not resolve'); } });
  fixture.state.similarityAlignmentPlan.value = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;
  fixture.state.linearRecordTranslations.value = [
    { recordKey: 'record-a', x: 1, y: 2 },
    { recordKey: 'record-b', x: 3, y: 4 }
  ];

  fixture.actions.retainForStableReorder(['record-b', 'record-a']);
  assert.deepEqual(
    fixture.state.similarityAlignmentPlan.value.records.map(({ recordKey }) => recordKey),
    ['record-b', 'record-a']
  );
  assert.deepEqual(
    fixture.state.linearRecordTranslations.value.map(({ recordKey }) => recordKey),
    ['record-b', 'record-a']
  );
});

test('source, crop, and selector invalidation publish their semantic reason', () => {
  const fixture = create({ helper: async () => { throw new Error('must not resolve'); } });
  const plan = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;
  for (const reason of ['source replaced.', 'record crop changed.', 'record selector changed.']) {
    fixture.state.similarityAlignmentPlan.value = structuredClone(plan);
    assert.equal(fixture.actions.clearForMutation(reason), true);
    assert.equal(fixture.state.similarityAlignmentPlan.value, null);
    assert.equal(fixture.actions.notice.value, `Alignment cleared: ${reason}`);
  }
});

test('ordinary Generate validation preserves a current plan without running generation or LOSATP', async () => {
  const fixture = create({
    helper: async (_operation, { request }) => ({ result: resolvedResponse(request) }),
    runAnalysis: async () => { throw new Error('validation must not generate'); }
  });
  fixture.state.similarityAlignmentPlan.value = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;

  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.state.similarityAlignmentPlan.value.groupId, 'og-1');
});

test('stale reference blocks Generate and offers explicit reference repair or Reset', async () => {
  const fixture = create({ helper: async () => { throw new Error('must not resolve'); } });
  const plan = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;
  plan.reference = anchor('record-a', 'removed-reference', 99);
  plan.records[0].anchor = plan.reference;
  // Vue exposes committed plans through a reactive Proxy; the controller's
  // canonical JSON boundary must not feed that Proxy to structuredClone.
  fixture.state.similarityAlignmentPlan.value = new Proxy(plan, {});

  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), {
    status: 'blocked', reason: 'stale-reference'
  });
  assert.equal(fixture.actions.repair.value.kind, 'reference');
  assert.match(fixture.actions.notice.value, /needs repair/);
});

test('stale target blocks Generate and requires Select or Skip without replacing Result', async () => {
  let call = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }) => {
      call += 1;
      if (call === 1) return { result: ambiguityResponse(request) };
      const records = [
        decision('record-a', 'reference', 'reference', request.reference),
        decision('record-b', 'skipped', 'skipped_by_user', null)
      ];
      return { result: {
        schema: 1, status: 'resolved', mode: request.mode, groupId: request.groupId,
        reference: request.reference, records,
        plan: {
          schema: 1, mode: request.mode, groupId: request.groupId,
          reference: request.reference, records: records.map(planDecision)
        }
      } };
    }
  });
  const plan = resolvedResponse({
    ...renderRequest,
    mode: 'position',
    groupId: 'og-1',
    reference: anchor('record-a', 'clicked-inparalog', 4)
  }).plan;
  plan.records[1].anchor = anchor('record-b', 'removed-target', 77);
  fixture.state.similarityAlignmentPlan.value = new Proxy(plan, {});

  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), {
    status: 'blocked', reason: 'stale-target'
  });
  assert.equal(fixture.actions.dialogOpen.value, true);
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
  assert.deepEqual(await fixture.actions.skipRecord('record-b'), { status: 'selected' });
});

test('failed validation and generation preserve local choices for retry', async () => {
  let validationAttempts = 0;
  let generationAttempts = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }, call) => {
      if (call === 1) return { result: ambiguityResponse(request) };
      validationAttempts += 1;
      if (validationAttempts === 1) throw new Error('Resolver temporarily unavailable');
      const response = resolvedResponse(request, request.choices[0].anchor);
      response.records[1].rationale = 'user_selected';
      response.plan.records[1].rationale = 'user_selected';
      return { result: response };
    },
    runAnalysis: async () => {
      generationAttempts += 1;
      return generationAttempts === 1
        ? { status: 'error', error: new Error('Generation temporarily unavailable') }
        : { status: 'ok' };
    }
  });
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  const selected = anchor('record-b', 'target-b', 2);
  fixture.actions.selectCandidate('record-b', selected);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'error' });
  assert.match(fixture.actions.error.value.message, /Resolver temporarily unavailable/);
  assert.equal(fixture.actions.canApply.value, true);
  assert.equal(generationAttempts, 0);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'error' });
  assert.match(fixture.actions.error.value.message, /Generation temporarily unavailable/);
  assert.equal(fixture.actions.canApply.value, true);
  assert.equal(fixture.actions.draft.value.ambiguities[0].choice.candidateKey, JSON.stringify(selected));
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 4);
  assert.equal(generationAttempts, 2);
  assert.equal(fixture.actions.draft.value, null);
});

test('changed committed Result rejects an old local draft before batch validation', async () => {
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }) => ({ result: ambiguityResponse(request) })
  });
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  fixture.actions.skipRecord('record-b');
  fixture.state.results.value = [{ name: 'newly-committed.svg' }];
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'stale' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.equal(fixture.actions.draft.value, null);
  assert.match(fixture.actions.error.value.message, /Start alignment again/);
});

test('candidate labels use biological metadata and displayed record names', async () => {
  const fixture = create({
    members: [reference, { ...targetA, gene: 'tnpA', locus_tag: 'ABC_00120' },
      { ...targetB, product: 'transposase' }],
    helper: async (_operation, { request }) => ({ result: ambiguityResponse(request) })
  });
  fixture.state.linearSeqs[1].definition = 'Bacillus chromosome B';
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  const record = fixture.actions.draft.value.ambiguities[0];
  assert.equal(record.recordLabel, 'Bacillus chromosome B');
  assert.equal(record.candidates[0].label, 'tnpA · ABC_00120');
  assert.equal(record.candidates[1].label, 'transposase');
  assert.equal(record.candidates[0].coordinates, '101..130 bp');
});

test('switching the selected Similarity Group invalidates an open draft', async () => {
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }) => ({ result: ambiguityResponse(request) })
  });
  fixture.state.selectedOrthogroupId = ref('og-1');
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  fixture.actions.skipRecord('record-b');
  fixture.state.selectedOrthogroupId.value = 'og-2';
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'stale' });
  assert.equal(fixture.helperCalls.length, 1);
});

test('Cancel during batch validation ignores the late Python response', async () => {
  let finishValidation;
  let generationCalls = 0;
  const fixture = create({
    members: [reference, targetA, targetB],
    helper: async (_operation, { request }, call) => call === 1
      ? { result: ambiguityResponse(request) }
      : new Promise((resolve) => { finishValidation = () => {
        const response = resolvedResponse(request, request.choices[0].anchor);
        response.records[1].rationale = 'user_selected';
        response.plan.records[1].rationale = 'user_selected';
        resolve({ result: response });
      }; }),
    runAnalysis: async () => { generationCalls += 1; return { status: 'ok' }; }
  });
  await fixture.actions.startFromPopup({ groupId: 'og-1', reference });
  fixture.actions.selectCandidate('record-b', anchor('record-b', 'target-a', 1));
  const applying = fixture.actions.applyDraft();
  while (!finishValidation) await Promise.resolve();
  fixture.actions.cancel();
  finishValidation();
  assert.deepEqual(await applying, { status: 'stale' });
  assert.equal(generationCalls, 0);
  assert.equal(fixture.actions.draft.value, null);
  assert.deepEqual(fixture.state.results.value, [{ name: 'committed.svg' }]);
});
