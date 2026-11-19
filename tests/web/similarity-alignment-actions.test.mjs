import assert from 'node:assert/strict';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import { drawingStore } from './helpers/drawing-state.mjs';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';

const ref = (value) => ({ value });
// The fixture refs are plain objects, so a watcher runs only on `flushWatchers()`.
const watchers = [];
const flushWatchers = () => watchers.forEach((entry) => {
  const value = entry.source();
  if (value !== entry.last) entry.callback((entry.last = value));
});
globalThis.window = {
  Vue: { ref, computed: (getter) => ({ get value() { return getter(); } }),
    watch: (source, callback) => { watchers.push({ source, callback, last: source() }); } }
};
const { createSimilarityAlignmentActions } = await import(
  '../../gbdraw/web/js/app/similarity-alignment.js'
);
const { createPasteTargets } = await import('../../gbdraw/web/js/app/region-drawing/paste-into.js');

const { buildSimilarityAlignmentResetReceipt } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
const { SOURCE_LOCALE, setUiLocale, t, tm } = await import('../../gbdraw/web/js/utils/i18n.js');
const { templateBinding } = await import('./helpers/template-binding.mjs');
const { errorText } = await import('../../gbdraw/web/js/utils/error-text.js');
const reviewView = await import('../../gbdraw/web/js/services/alignment-review-view.js');
const { writeCanonicalRecordReverseComplement } = await import('../../gbdraw/web/js/services/record-display-model.js');
const { CompositionMetadataError } = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');

const anchor = (recordKey, biologicalFeatureId, sourceFeatureIndex) => ({
  recordKey, biologicalFeatureId, sourceFeatureIndex,
  stableFeatureSvgId: `stable-${biologicalFeatureId}`
});
// The source facts of each fixture feature, as the committed request's records hold them: the fake helper reads
// them here by identity, as Python reads them from the request (a helper request names a member by identity only).
const sourceFacts = new Map();
const factsOf = ({ recordKey, biologicalFeatureId }) => sourceFacts.get(`${recordKey}\0${biologicalFeatureId}`);
const member = (recordKey, id, index, strand, representative = false) => {
  const start = 10 + index * 10;
  const end = 30 + index * 10;
  sourceFacts.set(`${recordKey}\0${id}`, { sourceStart: start, sourceEnd: end,
    sourceStrand: strand === '+' ? 1 : strand === '-' ? -1 : null });
  return {
    ...anchor(recordKey, id, index), recordIndex: recordKey === 'a' ? 0 : recordKey === 'b' ? 1 : 2,
    featureIndex: index, start, end,
    strand, proteinId: id, sourceProteinId: id, representative, role: 'member',
    gene: `gene-${id}`
  };
};
const reference = member('a', 'clicked', 4, '+');
const otherReference = member('a', 'other', 7, '+', true);
const targetA = member('b', 'target-a', 1, '-', true);
const targetB = member('b', 'target-b', 2, '+');
const targetC = member('c', 'target-c', 3, null);
const renderRequest = {
  schema: CANONICAL_REQUEST_SCHEMA, mode: 'linear', layout:{recordTranslations:[],similarityAlignment:null}, records: ['a', 'b', 'c'].map((recordKey) => ({
    recordKey, source:{kind:'genbank',resourceId:'source'}, selector:{kind:'recordId',id:recordKey},cardinality:'exactly_one',display:{isCircular:null,startCoordinate:null}, windows: [], presentation: { reverseComplement: false }
  }))
};
const state = () => {
  // The shown Linear drawing's rows.
  const linearSeqs = ['a', 'b', 'c'].map((uid) => ({ uid, region_reverse: false, windows: [] }));
  // What `sessionOperationAvailability('not-during-run')` returns: a test sets it while a run holds the owner.
  const runBlock = ref(null);
  return {
    runBlock,
    sessionOperationAvailability: (operation = 'mutation') => (operation === 'not-during-run' ? runBlock.value : null),
    ...drawingStore({ recordsOff: [] }, { recordsOff: [] }),
    activeDrawing: () => ({ mode: 'linear', recordsOff: [], linearSeqs }),
    activeDrawingId: ref('linear-1'),
    mode: ref('linear'),
    featureCatalog: ref({ items: [{ recordKeys: ['a', 'b', 'c'] }] }),
    selectedOrthogroupId: ref('og-1'),
    similarityAlignmentPlan: ref(null),
    similarityAlignmentResetReceipt:ref(null),
    errorLog: ref(null),
    linearRecordTranslations: ref([]),
    results: ref([{ name: 'prior.svg' }])
  };
};
const deferred = () => {
  let resolve;
  const promise = new Promise((done) => { resolve = done; });
  return { promise, resolve };
};
// A helper request record is drawn reversed by its presentation or by its windows, all reversed.
const reversedRow = (row) => row.presentation.reverseComplement
  || (row.windows.length > 0 && row.windows.every((window) => window.reverseComplement));
const serializedCandidate = (member, request) => {
  const item = { ...member, ...factsOf(member.anchor) };
  const reverse = key => reversedRow(request.records.find(r=>r.recordKey===key));
  const displayedStrand = item.sourceStrand === null ? null : item.sourceStrand*(reverse(item.anchor.recordKey)?-1:1);
  const referenceSourceStrand = factsOf(request.reference)?.sourceStrand;
  const referenceStrand=referenceSourceStrand===null?null:referenceSourceStrand*(reverse(request.reference.recordKey)?-1:1);
  const strandRelation = displayedStrand === null || referenceStrand === null
    ? 'unknown' : displayedStrand === referenceStrand ? 'same' : 'opposite';
  return {
    anchor: item.anchor, displayName: item.anchor.biologicalFeatureId,
    sourceStart: item.sourceStart, sourceEnd: item.sourceEnd,
    displayCenter: reverse(item.anchor.recordKey)?300-(item.sourceStart + item.sourceEnd)/2:(item.sourceStart + item.sourceEnd)/2,
    displayedStrand, hidden: false, representative: item.representative,
    role: item.role, usable: true, directEvidence: [], strandRelation
  };
};
const responseFor = (request, { ambiguous = false, missing = false } = {}) => {
  const records = request.records.map(({ recordKey }) => {
    if (recordKey === request.reference.recordKey) return {
      kind: 'decision', recordKey, status: 'reference', rationale: 'reference',
      reviewReason: null, anchor: request.reference,
      candidates: []
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
        candidates
      };
    }
    if (!candidates.length || missing && recordKey === 'c') return {
      kind: 'decision', recordKey, status: 'skipped',
      rationale: 'skipped_no_candidate', reviewReason: null, anchor: null,
      candidates: []
    };
    if (ambiguous && recordKey === 'b') return {
      kind: 'ambiguous', recordKey, candidates, directRbhCandidates: [],
      recommendedAnchor: candidates[0].anchor,
      recommendationReason: 'unique_representative'
    };
    return {
      kind: 'decision', recordKey, status: 'aligned',
      rationale: 'only_usable_candidate', reviewReason: 'only_usable_candidate',
      anchor: candidates[0].anchor, candidates
    };
  });
  const status = records.some(({ kind }) => kind === 'ambiguous') ? 'ambiguous' : 'resolved';
  return {
    schema: 2, status, groupId: request.groupId, reference: request.reference,
    referenceDisplayedStrand: 1, referenceDisplayCenter: 65,
    referenceSourceStart: factsOf(request.reference).sourceStart, referenceSourceEnd: factsOf(request.reference).sourceEnd,
    records, projection: null,
    plan: status === 'resolved' ? {
      schema: 2, groupId: request.groupId, reference: request.reference,
      records: records.map(({ recordKey, status: decisionStatus, rationale, anchor: chosen }) => ({
        recordKey, status: decisionStatus, rationale, anchor: chosen
      }))
    } : null
  };
};
const create = ({
  members = [reference, targetA, targetC], helper = (_operation, { request }) => ({
    result: responseFor(request)
  }), generation = async () => ({ status: 'ok' }), onError = null,
  previewCandidate = null, clearCandidatePreview = null, mutateProjection = null,
  getCommittedRequest = null, writeArtifactValues = () => {}, drawingNameOf = (drawingId) => `Drawing ${drawingId}`
} = {}) => {
  const controllerState = state();
  controllerState.featureCatalog.value.items[0].biologicalFeatures=members;
  const currentGroup = { id: 'og-1', members, orthologEdges: [] };
  let enrichCalls = 0;
  const helperCalls = [];
  const generationCalls = [];
  let svg = null;
  const bytes=Buffer.from('fixture source');
  let canonical={renderRequest:structuredClone(renderRequest),resources:{source:{name:'source.gbk',kind:'genbank',encoding:'base64',data:bytes.toString('base64'),size:bytes.length}}};
  // Counts replacements of the committed Session, as services/config.js does, so a computed follows them.
  const committedRevision = ref(0);
  const projectCommittedAlignment=({committed,plan,translations,orientations})=>{
    const next={...committed,renderRequest:structuredClone(committed.renderRequest)};
    next.renderRequest.layout={recordTranslations:structuredClone(translations),similarityAlignment:structuredClone(plan??null)};
    for(const direction of orientations)writeCanonicalRecordReverseComplement(next.renderRequest.records.find(r=>r.recordKey===direction.recordKey),direction.reverseComplement);
    return next;
  };
  // History compares handles by their identity (OV-2208): this fake's
  // fingerprint is its owner set.
  const capture = () => {
    const ownerSet = structuredClone({
      canonical, receipt:controllerState.similarityAlignmentResetReceipt.value,
      plan: controllerState.similarityAlignmentPlan.value,
      translations: controllerState.linearRecordTranslations.value,
      sequences: controllerState.activeDrawing().linearSeqs,
      results: controllerState.results.value
    });
    return { retainedBytes: 0, identity: { fingerprint: JSON.stringify(ownerSet), compactSignature: '' }, ownerSet };
  };
  const history = createHistoryManager({
    buildIntent: () => ({}), applyIntent: () => {},
    buildCheckpoint: () => { throw new Error('Apply must use generated artifact handles.'); },
    applyCheckpoint: () => {},
    captureGeneratedArtifactHandle: capture,
    restoreGeneratedArtifactHandle: ({ ownerSet }) => {
      canonical=structuredClone(ownerSet.canonical);
      controllerState.similarityAlignmentResetReceipt.value=structuredClone(ownerSet.receipt);
      controllerState.similarityAlignmentPlan.value = structuredClone(ownerSet.plan);
      controllerState.linearRecordTranslations.value = structuredClone(ownerSet.translations);
      controllerState.activeDrawing().linearSeqs.splice(0, controllerState.activeDrawing().linearSeqs.length,
        ...structuredClone(ownerSet.sequences));
      controllerState.results.value = structuredClone(ownerSet.results);
    }
  });
  const actions = createSimilarityAlignmentActions({
    state: controllerState,
    writeArtifactValues,
    drawingNameOf,
    getOrthogroupById: (id) => id === currentGroup.id ? currentGroup : null,
    getEnrichedOrthogroupMembers: () => { enrichCalls += 1; return currentGroup.members; },
    // The orthogroup owner's drawn IDs of an anchor's feature; this Result names its records by key.
    drawnFeatureIdsOf: () => [],
    // The orthogroup owner's detail of the selected group (OV-920).
    selectedOrthogroupDetail: { get value() {
      return { group: currentGroup, id: currentGroup.id, members: currentGroup.members };
    } },
    getCommittedRequest:()=>{ void committedRevision.value; return getCommittedRequest?.()||canonical.renderRequest; },
    getCommittedSession:()=>({...canonical,renderRequest:getCommittedRequest?.()||canonical.renderRequest}),
    projectCommittedAlignment,
    getCurrentSvg: () => svg,
    // The composition root's port: run the candidate, then commit the changed orientations.
    runRecordAlignment: (options) => history.runUndoableArtifactReplacement(options.label, async () => {
      generationCalls.push(options);
      const outcome = await generation(options);
      if (outcome.error?.summary) controllerState.errorLog.value = outcome.error;
      if (outcome.status === 'ok') {
        const receipt=options.alignmentResetBefore?await buildSimilarityAlignmentResetReceipt({before:options.alignmentResetBefore,after:options.canonical}):options.alignmentResetReceipt;
        canonical=options.canonical;
        controllerState.similarityAlignmentResetReceipt.value=receipt;
        controllerState.similarityAlignmentPlan.value = canonical.renderRequest.layout.similarityAlignment;
        controllerState.linearRecordTranslations.value = canonical.renderRequest.layout.recordTranslations;
        options.orientations.forEach(({recordKey,reverseComplement})=>{controllerState.activeDrawing().linearSeqs.find(r=>r.uid===recordKey).region_reverse=reverseComplement;});
        controllerState.results.value = [{ name: 'aligned.svg' }];
        svg = null;
      }
      return outcome;
    }, { shouldCommit: (outcome) => outcome.status === 'ok' }),
    cancelRunAnalysis: () => {},
    runHelperOperation: async (operation, payload) => {
      helperCalls.push({ operation, payload: structuredClone(payload) });
      const output=await helper(operation, payload, helperCalls.length);
      if(output?.result && payload.projection && payload.projection.geometry!==false && output.result.projection===null){
        const r=output.result;const request=payload.request;
        const current=payload.projection.canonicalRequest;
        const refMember=request.members.find(m=>JSON.stringify(m.anchor)===JSON.stringify(request.reference));
        const refCandidate=serializedCandidate(refMember,request);
        r.referenceDisplayCenter=refCandidate.displayCenter;
        r.referenceDisplayedStrand=refCandidate.displayedStrand;
        r.projection={binding:'a'.repeat(64),geometryOrientations:payload.projection.orientations||Object.fromEntries(request.records.map(r=>[r.recordKey,Boolean(r.presentation?.reverseComplement)])),records:request.records.map(record=>{
          const candidates=request.members.filter(m=>m.anchor.recordKey===record.recordKey).map(m=>({...m,...factsOf(m.anchor)}));
          const before=reversedRow(record);
          const facts=reverse=>candidates.map(m=>{
            const candidate=r.records.find(row=>row.recordKey===record.recordKey).candidates.find(c=>JSON.stringify(c.anchor)===JSON.stringify(m.anchor));
            const center=candidate?.usable===false?null:reverse?300-(m.sourceStart+m.sourceEnd)/2:(m.sourceStart+m.sourceEnd)/2;
            return {anchor:m.anchor,displayCenter:center,centerX:center,displayedStrand:m.sourceStrand===null?null:m.sourceStrand*(reverse?-1:1)};
          });
          return {recordKey:record.recordKey,beforeReverseComplement:before,base:current.layout?.recordTranslations.find(t=>t.recordKey===record.recordKey)||{x:0,y:0},beforeAxisY:100,
            beforeAnchors:facts(before).filter(f=>f.centerX!==null).map(({anchor,centerX})=>({anchor,centerX})),variants:[false,true].map(reverseComplement=>({reverseComplement,axisY:100,anchors:facts(reverseComplement)}))};
        })};
        r.projection.records.forEach(f=>{delete f.base.recordKey;});
      }
      mutateProjection?.(output.result,helperCalls.length);
      return output;
    },
    resolveOperation: 'resolveSimilarityAlignment',
    onError: (error) => { controllerState.errorLog.value = error; onError?.(error); },
    previewCandidate, clearCandidatePreview
  });
  return { actions, state: controllerState, helperCalls, generationCalls, currentGroup, history,
    enrichCalls: () => enrichCalls,
    // A Session load or a Generate of the same files commits a Session with its own resources.
    replaceCommittedResources: () => {
      canonical = { ...canonical, resources: structuredClone(canonical.resources) };
      committedRevision.value += 1;
    },
    showResultSvg: (value) => { svg = value; } };
};
const startReview = (fixture, selected = reference) => fixture.actions.startFromPopup({
  groupId: 'og-1', reference: selected, mode: 'review'
});
const startAlign = (fixture, selected = reference) => fixture.actions.startFromPopup({
  groupId: 'og-1', reference: selected
});

// OV-2771: Python places an anchor by the window that holds it, so the helper request carries each record's
// windows as the request stores them, a window across the origin included.
test('the helper request carries each record\'s windows as stored, a window across the origin included', async () => {
  const committed = structuredClone(renderRequest);
  committed.records[1].windows = [{ start: 250, end: 40, reverseComplement: true, column: null }];
  committed.records[2].windows = [{ start: 11, end: 200, reverseComplement: false, column: null }];
  const fixture = create({ members: [reference, targetA, targetC], getCommittedRequest: () => committed });
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
  const [{ payload }] = fixture.helperCalls;
  assert.deepEqual(payload.request.records.map(({ recordKey, windows, presentation }) => ({ recordKey, windows, presentation })), [
    { recordKey: 'a', windows: [], presentation: { reverseComplement: false } },
    { recordKey: 'b', windows: [{ start: 250, end: 40, reverseComplement: true, column: null }],
      presentation: { reverseComplement: false } },
    { recordKey: 'c', windows: [{ start: 11, end: 200, reverseComplement: false, column: null }],
      presentation: { reverseComplement: false } }
  ]);
  assert.deepEqual(payload.projection.orientations, { a: false, b: true, c: false });
});

// D5: Python reads each record's length from the committed request, so the page sends none.
test('the helper request record carries no record length', async () => {
  const fixture = create();
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
  const [{ payload }] = fixture.helperCalls;
  assert.deepEqual(payload.request.records.map((record) => Object.keys(record).sort()),
    [0, 1, 2].map(() => ['presentation', 'recordKey', 'windows']));
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
      unusable.candidates[0].usable = false;
      unusable.candidates[0].displayCenter = null;
      Object.assign(response.plan.records[2], {
        status: 'skipped', rationale: 'skipped_unmappable', anchor: null
      });
      return { result: response };
    },
    generation: () => applying.promise });
  const priorResult = fixture.state.results.value;
  const operation = startAlign(fixture);
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(fixture.actions.status.value, 'applying');
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.ok(fixture.actions.draft.value);
  assert.equal(fixture.actions.summary.value, null);
  assert.equal(fixture.state.results.value, priorResult);
  assert.equal(fixture.helperCalls.length, 1);
  assert.deepEqual(await startAlign(fixture), { status: 'busy', reason: 'Applying the alignment plan.' });
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

// OV-1010: a member's null or blank bound was read as 0. The page now sends no bound: Python reads the member's
// source coordinates from the committed request, so a saved member without them aligns where its feature is.
for (const missing of [null, '', '  ']) {
  test(`Align reads a member whose saved start is ${JSON.stringify(missing)} from Python (OV-1010)`, async () => {
    const fixture = create({ members: [reference, { ...targetA, start: missing }, targetC] });
    assert.deepEqual(await startAlign(fixture), { status: 'ok' });
    assert.equal(fixture.helperCalls.length, 1);
    assert.equal('sourceStart' in fixture.helperCalls[0].payload.request.members[1], false);
  });
}

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
      // The Worker reads direct evidence from the committed resource; a catalog graph is ignored.
      assert.deepEqual(request.directEdges, []);
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
  assert.equal(fixture.actions.draft.value.rows[0].choice, null);
  assert.equal(fixture.actions.canApply.value, false);
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
  const renderError = { code: 'COMPARISON_IDENTITY', stage: 'render', context: { reason: 'SOURCE_VIEW_CONFLICT' }, message: 'PRIVATE_COMPARE_SENTINEL' };
  const fixture = create({ generation: async () => fail
    ? { status: 'error', error: renderError }
    : { status: 'ok' } });
  const previous = fixture.state.results.value;
  assert.equal((await startAlign(fixture)).status, 'error');
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.dialogOpen.value, true);
  assert.equal(fixture.actions.draft.value.response.status, 'resolved');
  assert.equal(fixture.actions.error.value.code, 'COMPARISON_IDENTITY');
  assert.equal(fixture.actions.error.value.operation, 'align');
  assert.equal(fixture.actions.error.value.context.reason, 'SOURCE_VIEW_CONFLICT');
  assert.doesNotMatch(JSON.stringify(fixture.actions.error.value), /PRIVATE_COMPARE_SENTINEL/);
  assert.equal(fixture.state.errorLog.value, fixture.actions.error.value);
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

// While a Generate runs, the run refuses an Align and a Reset (`sessionOperationAvailability('not-during-run')`):
// a refusal, not a failure, so no error is published and the Result stays (OV-2446).
test('an Align or Reset run refused while a Generate runs publishes no error and keeps the Result', async () => {
  const busy = { status: 'busy', reason: 'Generating diagram. Retry after generation finishes.' };
  let refuse = true;
  const f = create({ generation: async () => (refuse ? busy : { status: 'ok' }) });
  const prior = f.state.results.value;
  assert.deepEqual(await startAlign(f), busy);
  assert.equal(f.state.errorLog.value, null);
  assert.equal(f.actions.error.value, null);
  assert.equal(f.state.results.value, prior);
  assert.equal(f.history.getUndoCount(), 0);
  f.actions.cancel();
  refuse = false;
  await startReview(f); f.actions.setDirectionMode('left'); assert.equal((await f.actions.applyDraft()).status, 'ok');
  const aligned = f.state.results.value;
  refuse = true;
  assert.deepEqual(await f.actions.resetAlignment('positions'), busy);
  assert.equal(f.state.errorLog.value, null);
  assert.equal(f.actions.error.value, null);
  assert.equal(f.state.results.value, aligned);
  assert.equal(f.history.getUndoCount(), 1);
});

test('busy is synchronous, duplicate starts are ignored, and exact popup reference survives', async () => {
  const pending = deferred();
  const fixture = create({ members: [reference, otherReference, targetA, targetC],
    helper: (_operation, { request }) => pending.promise.then(() => ({ result: responseFor(request) })) });
  const operation = startReview(fixture);
  assert.equal(fixture.actions.status.value, 'resolving');
  assert.equal(fixture.actions.busy.value, true);
  assert.deepEqual(await startReview(fixture, otherReference), { status: 'busy', reason: 'Resolving the current alignment.' });
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
    assert.equal(fixture.actions.canApply.value, !ambiguous);
    assert.equal(fixture.actions.unresolvedCount.value, Number(ambiguous));
    if(ambiguous)fixture.actions.selectCandidate('b',anchor('b','target-a',1));
    assert.deepEqual(fixture.actions.draft.value.rows.map(({ recordKey }) => recordKey), ['b', 'c']);
    assert.equal(fixture.actions.draft.value.rows[0].choice.candidateKey,
      JSON.stringify(anchor('b', 'target-a', 1)));
    assert.equal(fixture.actions.draft.value.rows[0].reason,
      ambiguous ? 'user_selected' : 'only_usable_candidate');
    // Names come from biological metadata; the feature hash stays binding-only.
    assert.equal(reviewView.alignmentFeatureName(fixture.actions.draft.value.reference.feature).label, 'gene-clicked · clicked');
    assert.equal(Object.hasOwn(fixture.actions.draft.value.reference, 'featureIdentifier'), false);
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
    [recordKey, kind]), [['b', 'select']]);
  assert.equal(fixture.state.similarityAlignmentPlan.value.records.find(row => row.recordKey === 'c').rationale, 'skipped_no_candidate');
  assert.equal(fixture.state.similarityAlignmentPlan.value.schema, 2);
  assert.deepEqual(fixture.state.results.value, [{ name: 'aligned.svg' }]);
});

test('candidate selection and Skip are local edits with zero Worker jobs', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  await startReview(fixture);
  const count = fixture.helperCalls.length;
  assert.deepEqual(fixture.actions.selectCandidate('b', anchor('b', 'target-b', 2)),
    { status: 'selected' });
  assert.equal(fixture.actions.draft.value.rows[0].choice.candidateKey,
    JSON.stringify(anchor('b', 'target-b', 2)));
  assert.equal(fixture.actions.draft.value.rows[0].candidates[1].strandRelation, 'same');
  fixture.actions.skipRecord('b');
  assert.equal(fixture.actions.draft.value.rows[0].choice.kind, 'skip');
  fixture.actions.selectCandidate('b', anchor('b', 'target-a', 1));
  assert.equal(fixture.actions.draft.value.rows[0].candidates[0].strandRelation, 'opposite');
  assert.equal(fixture.actions.draft.value.rows[1].candidates[0].strandRelation, 'unknown');
  assert.equal(fixture.helperCalls.length, count);
  assert.equal(fixture.generationCalls.length, 0);
});

test('Apply validates one explicit batch and commits through one generation call', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  await startReview(fixture);
  fixture.actions.selectCandidate('b',anchor('b','target-a',1));
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 1);
  assert.deepEqual(fixture.helperCalls[1].payload.request.choices, [
    { recordKey: 'b', kind: 'select', anchor: anchor('b', 'target-a', 1) },
    { recordKey: 'c', kind: 'select', anchor: anchor('c', 'target-c', 3) }
  ]);
  assert.equal(fixture.actions.summary.value.reversed, 0);
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
  // The work in flight refuses both, with its reason.
  const applyingWork = { status: 'busy', reason: 'Applying the alignment plan.' };
  assert.deepEqual(await fixture.actions.applyDraft(), applyingWork);
  assert.deepEqual(await startReview(fixture), applyingWork);
  pending.resolve();
  assert.deepEqual(await operation, { status: 'ok' });
  assert.equal(fixture.helperCalls.length, 2);
  assert.equal(fixture.generationCalls.length, 1);
});

test('validation and generation failures preserve editable draft and prior artifact', async () => {
  for (const failure of ['validation', 'generation']) {
    let fail = true;
    const renderError = { code: 'COMPARISON_IDENTITY', stage: 'render', context: { reason: 'SOURCE_VIEW_CONFLICT' }, message: 'PRIVATE_COMPARE_SENTINEL' };
    const fixture = create({
      helper: (_operation, { request }, count) => {
        if (failure === 'validation' && count > 1 && fail) throw new Error('Invalid choice. Select another.');
        return { result: responseFor(request) };
      },
      generation: async () => failure === 'generation' && fail
        ? { status: 'error', error: renderError }
        : { status: 'ok' }
    });
    await startReview(fixture);
    fixture.actions.setDirectionMode('right');
    const before = fixture.actions.draft.value;
    const orientations = fixture.state.activeDrawing().linearSeqs.map(({ region_reverse }) => region_reverse);
    assert.equal((await fixture.actions.applyDraft()).status, 'error');
    assert.deepEqual(fixture.state.activeDrawing().linearSeqs.map(({ region_reverse }) => region_reverse), orientations);
    assert.equal(fixture.history.getUndoCount(), 0);
    assert.equal(fixture.actions.status.value, 'reviewing');
    assert.equal(fixture.actions.draft.value, before);
    assert.deepEqual(fixture.state.results.value, [{ name: 'prior.svg' }]);
    assert.equal(fixture.state.similarityAlignmentPlan.value, null);
    if (failure === 'generation') {
      assert.equal(fixture.actions.error.value.code, 'COMPARISON_IDENTITY');
  assert.equal(fixture.actions.error.value.operation, 'align');
  assert.equal(fixture.actions.error.value.context.reason, 'SOURCE_VIEW_CONFLICT');
  assert.doesNotMatch(JSON.stringify(fixture.actions.error.value), /PRIVATE_COMPARE_SENTINEL/);
      assert.equal(fixture.state.errorLog.value, fixture.actions.error.value);
    } else {
      assert.equal(fixture.actions.error.value.code, 'UNKNOWN');
      assert.doesNotMatch(fixture.actions.error.value.summary, /Invalid choice/);
    }
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

// Delta review B1 of #2047: Copy to another drawing into the drawing of an open review adds records the review did
// not read; its Apply would draw the committed records only and leave the drawing's plan short of the pasted ones.
// The paste's History step ends that drawing's alignment work through `endWorkOf`, as a record edit does, and keeps
// the plan and translations the paste installs (Python's); a review of another drawing stays (D-38).
test('a paste into the drawing of an open review ends the review and keeps the plan the paste installs', async () => {
  for (const aligned of [false, true]) {
    const fixture = create();
    const own = { schema: 2, groupId: 'og-1', reference: anchor('a', 'clicked', 4),
      records: [{ recordKey: 'a', status: 'reference', rationale: 'reference', anchor: anchor('a', 'clicked', 4) }] };
    if (aligned) {
      fixture.state.similarityAlignmentPlan.value = structuredClone(own);
      fixture.state.linearRecordTranslations.value = [{ recordKey: 'a', x: 30, y: 0 }];
    }
    await startReview(fixture);
    assert.equal(fixture.actions.status.value, 'reviewing');
    const installed = aligned
      ? { similarityAlignmentPlan: { ...own, records: [...own.records,
        { recordKey: 'pasted', status: 'skipped', rationale: 'skipped_no_candidate', anchor: null }] },
      linearRecordTranslations: [{ recordKey: 'a', x: 30, y: 0 }, { recordKey: 'pasted', x: 0, y: 0 }] }
      : { similarityAlignmentPlan: null, linearRecordTranslations: [] };
    const pasteInto = (drawingId) => createPasteTargets({
      tabs: () => [], drawingOf: (id) => ({ id }), readRecords: async () => undefined, project: async () => null,
      show: () => ({ status: 'ok' }), runUndoable: async (_label, run) => run(), prepare: async () => ({}),
      install: () => structuredClone(installed),
      // Only the review's drawing (linear-1) is this fixture's state; a write for another drawing goes elsewhere.
      writeArtifactValues: (writtenId, values) => Object.entries(values).forEach(([key, value]) => {
        if (writtenId === 'linear-1') fixture.state[key].value = value;
      }),
      endAlignmentWorkOf: fixture.actions.endWorkOf
    }).pasteInto(drawingId, {});
    const notice = fixture.actions.notice.value;
    assert.equal((await pasteInto('linear-2')).status, 'ok');
    assert.equal(fixture.actions.status.value, 'reviewing', 'a paste into another drawing keeps this review');
    if (aligned) {
      assert.notDeepEqual(fixture.state.similarityAlignmentPlan.value, installed.similarityAlignmentPlan,
        'a paste into another drawing writes nothing here');
    }
    assert.equal((await pasteInto('linear-1')).status, 'ok');
    assert.equal(fixture.actions.status.value, 'idle', 'the review of the pasted drawing ends');
    assert.equal(fixture.actions.draft.value, null);
    assert.equal(fixture.actions.dialogOpen.value, false);
    assert.equal(fixture.actions.canApply.value, false);
    assert.deepEqual(await fixture.actions.applyDraft(), { status: 'rejected' }, 'no Apply of the records it read');
    assert.equal(fixture.generationCalls.length, 0);
    assert.deepEqual(fixture.state.similarityAlignmentPlan.value, installed.similarityAlignmentPlan, 'the plan is kept');
    assert.deepEqual(fixture.state.linearRecordTranslations.value, installed.linearRecordTranslations);
    assert.equal(fixture.actions.notice.value, notice, 'no plan-cleared notice');
  }
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

test('a new Result or catalog under the same committed resources keeps the draft; other resources reject it before validation', async () => {
  const fixture = create();
  await startReview(fixture);
  // A live setting change draws a new Result from the same sources (PD-OI-102, OV-1082), and
  // the automatic rerender a new feature catalog of the same features (OV-1120). It keeps the
  // committed resources.
  fixture.state.results.value = [{ name: 'new.svg' }];
  const catalog = fixture.state.featureCatalog.value;
  fixture.state.featureCatalog.value = { ...catalog, items: catalog.items.map((item) => ({ ...item })) };
  assert.equal(fixture.actions.isDraftArtifactCurrent(), true);
  // A Generate or a Session load, even of the same files and features, commits other resources:
  // it replaces the figure (PD-OI-039). The resources decide, not the catalog.
  fixture.replaceCommittedResources();
  assert.equal(fixture.actions.isDraftArtifactCurrent(), false);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'stale' });
  assert.equal(fixture.helperCalls.length, 1);
  assert.deepEqual(fixture.state.results.value, [{ name: 'new.svg' }]);
});

// OV-1120: the automatic rerender commits a live setting during the review. The helper, Apply
// and the Reset receipt read the committed Session as it is then (PD-OI-102).
test('Apply composes onto the committed Session that a live setting changed during the review', async () => {
  let committed = structuredClone(renderRequest);
  const fixture = create({ getCommittedRequest: () => committed });
  await startReview(fixture);
  committed = { ...structuredClone(committed), diagramOptions: { plotTitle: 'live' } };
  assert.equal(fixture.actions.isDraftArtifactCurrent(), true);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.helperCalls.at(-1).payload.projection.canonicalRequest.diagramOptions.plotTitle, 'live');
  const [call] = fixture.generationCalls;
  assert.equal(call.canonical.renderRequest.diagramOptions.plotTitle, 'live');
  assert.equal(call.alignmentResetBefore.renderRequest.diagramOptions.plotTitle, 'live');
});

test('active schema-2 plan validates on regeneration and supports Reset and stable reorder', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.activePlanInspector.value.summary, 'Aligned to gene-clicked in Record 1');
  const calls = fixture.helperCalls.length;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, calls + 1);
  // The check resolves on the committed request, as every helper call does (D5-13: the projection is required), and
  // asks for the decisions only: Python plans and builds no layout for it.
  const { payload: check } = fixture.helperCalls.at(-1);
  assert.equal(check.projection.canonicalRequest.records.length, 3);
  assert.equal(check.projection.orientations, null);
  assert.equal(check.projection.geometry, false);
  assert.ok(check.resources.source);
  fixture.actions.retainForStableReorder(['c', 'a', 'b']);
  assert.deepEqual(fixture.state.similarityAlignmentPlan.value.records.map(({ recordKey }) => recordKey),
    ['c', 'a', 'b']);
  assert.deepEqual((await fixture.actions.resetAlignment()).status, 'ok');
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
});

// The check before Generate asks Python for decisions only, so it reads no record placement off the shown Result: a
// Result an older gbdraw drew (a Session 27 to 44) has no record keys to read placements by.
test('the check before Generate reads nothing from the shown Result', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.showResultSvg(new Proxy({}, { get: () => { throw new Error('The check read the shown Result.'); } }));
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  const { payload: check } = fixture.helperCalls.at(-1);
  assert.equal(check.projection.geometry, false);
  assert.equal(check.projection.canonicalRequest.layout.similarityAlignment, null);
});

// A Session an older gbdraw wrote saved a member's strand in the drawn direction: the page names a member by its
// identity, and Python reads its source facts from the committed request.
test('a helper request names each group member by its identity only', async () => {
  const fixture = create({ members: [reference, targetA, targetC].map((item) => ({ ...item })) });
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.currentGroup.members.forEach((item) => { item.strand = item.strand === '-' ? '+' : item.strand; });
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  for (const { payload } of fixture.helperCalls) {
    for (const item of payload.request.members) {
      assert.deepEqual(Object.keys(item).sort(),
        ['anchor', 'groupId', 'hidden', 'identityIsUnique', 'representative', 'role']);
    }
  }
});

test('Python recommendation remains a hint while direct-RBH resolved selection is preserved', async () => {
  const ambiguous = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => {
      const response = responseFor(request, { ambiguous: true });
      response.records[1].recommendedAnchor = anchor('b', 'target-b', 2);
      response.records[1].recommendationReason = 'deterministic_candidate_1';
      return { result: response };
    } });
  await startReview(ambiguous);
  assert.equal(ambiguous.actions.draft.value.rows[0].choice,null);
  assert.equal(ambiguous.actions.draft.value.rows[0].recommendedKey,
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

// OV-3440: a Result an older gbdraw drew has no record keys to read a stale plan's placements by. The repair of the
// check before Generate then places the records where the committed request has them, as the next Generate draws
// them, instead of stopping Generate.
test('the stale-target repair starts from the committed placements when the shown Result cannot be read', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const committedTranslations = structuredClone(fixture.state.linearRecordTranslations.value);
  fixture.showResultSvg(new Proxy({}, { get: () => {
    throw new CompositionMetadataError('Composition metadata cannot bind alignment record "b".');
  } }));
  fixture.currentGroup.members.splice(1, 1, targetB);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-target' });
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.repair.value.kind, 'targets');
  const { payload: repair } = fixture.helperCalls.at(-1);
  assert.deepEqual(repair.projection.canonicalRequest.layout.recordTranslations, committedTranslations);
  assert.equal(repair.projection.canonicalRequest.layout.similarityAlignment, null);
  assert.equal(fixture.actions.selectCandidate('b', anchor('b', 'target-b', 2)).status, 'selected');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.state.similarityAlignmentPlan.value.records[1].anchor.biologicalFeatureId, 'target-b');
});

// OV-3538: Reset record positions replaces every placement, so a Result an older gbdraw drew, which names no record by
// key, gives it the committed Session through the stale-target repair's rule instead of refusing it.
test('the repair base Session reads the committed placements when the shown Result cannot be read', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const committedTranslations = structuredClone(fixture.state.linearRecordTranslations.value);
  fixture.showResultSvg(new Proxy({}, { get: () => {
    throw new CompositionMetadataError('Composition metadata cannot bind alignment record "b".');
  } }));
  assert.throws(() => fixture.actions.materializedSession(), CompositionMetadataError);
  const { renderRequest } = fixture.actions.repairBaseSession();
  assert.deepEqual(renderRequest.layout.recordTranslations, committedTranslations);
  assert.equal(renderRequest.layout.similarityAlignment, null);
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

// OV-920: the drawer of the selected group reads the members the orthogroup
// owner enriched for it; reading its options and reason enriches none again.
test('the selected group drawer reads its members from the selected detail (OV-920)', () => {
  const fixture = create({ members: [reference, otherReference, targetA, targetC] });
  const before = fixture.enrichCalls();
  const drawer = fixture.actions.drawer.value;
  assert.equal(drawer.groupId, 'og-1');
  assert.deepEqual(drawer.options.map(({ label }) => label), fixture.actions.drawerReferenceOptions('og-1').map(({ label }) => label));
  assert.equal(drawer.options.length, 4);
  assert.match(drawer.disabledReason, /Select an exact reference/);
  assert.equal(fixture.actions.drawerDisabledReason('og-1'), drawer.disabledReason);
  assert.equal(fixture.enrichCalls(), before, 'the selected group is not enriched again');
  assert.equal(fixture.actions.setDrawerReference('og-1', drawer.options[1].key), true);
  assert.equal(fixture.actions.drawer.value.disabledReason, '');
  assert.deepEqual(fixture.actions.drawerReferenceOptions('og-missing'), []);
});

// OV-3442: a Session an older gbdraw wrote saved a member's strand in the drawn direction (a record drawn
// reverse-complemented). Two members with one name in one record show their location, read from the source facts their
// name reads: the catalog row after a Generate, the member's catalog feature (`sequenceFeature`) before one.
test('the drawer names an inparalog by the location and strand of its source feature', () => {
  const inparalogs = [member('b', 'para-1', 5, '-'), member('b', 'para-2', 6, '-')]
    .map((item) => ({ ...item, gene: 'gene-para' }));
  const fixture = create({ members: [reference, ...inparalogs] });
  // The group's members as such a Session saved them; the catalog keeps the source facts.
  fixture.currentGroup.members = fixture.currentGroup.members.map((item) => (
    item.strand === '-' ? { ...item, strand: '+', start: item.start + 1, sequenceFeature: item } : item
  ));
  const labels = () => fixture.actions.drawer.value.options.map(({ label }) => label).slice(1);
  const expected = ['Record 2 · gene-para · 61..80 (-)', 'Record 2 · gene-para · 71..90 (-)'];
  assert.deepEqual(labels(), expected);
  // A loaded Session before its first Generate has no catalog rows yet.
  fixture.state.featureCatalog.value = null;
  assert.deepEqual(labels(), expected);
});

// OV-381: a committed plan (from the popup or a loaded Session) shows its exact
// reference in the drawer; the drawer offers that reference for Align until the
// user picks another one, and never offers a stale one.
test('drawer offers the active plan exact reference until the user picks another (OV-381)', async () => {
  const fixture = create({ members: [reference, otherReference, targetA, targetC] });
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const planKey = JSON.stringify(fixture.state.similarityAlignmentPlan.value.reference);
  assert.equal(fixture.actions.drawerReferenceKey.value, planKey);
  assert.equal(fixture.actions.drawerDisabledReason('og-1'), '');
  const calls = fixture.helperCalls.length;
  assert.deepEqual(await fixture.actions.startFromDrawer({ groupId: 'og-1', mode: 'review' }),
    { status: 'reviewing' });
  assert.deepEqual(fixture.helperCalls[calls].payload.request.reference,
    fixture.state.similarityAlignmentPlan.value.reference);
  fixture.actions.cancel();
  const other = fixture.actions.drawerReferenceOptions('og-1')
    .find(({ anchor: a }) => a.biologicalFeatureId === 'other');
  assert.equal(fixture.actions.setDrawerReference('og-1', other.key), true);
  assert.equal(fixture.actions.drawerReferenceKey.value, other.key);
  // Choosing "none" explicitly is kept as well.
  assert.equal(fixture.actions.setDrawerReference('og-1', ''), false);
  assert.match(fixture.actions.drawerDisabledReason('og-1'), /Select an exact reference/);
});

test('drawer does not offer a stale plan reference (OV-381)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.currentGroup.members.splice(0, 1);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.actions.drawerReferenceKey.value, '');
  assert.notEqual(fixture.actions.drawerDisabledReason('og-1'), '');
});

test('a new committed plan drops the drawer pick (OV-381)', async () => {
  const planKey = (fixture) => JSON.stringify(fixture.state.similarityAlignmentPlan.value.reference);
  // Drawer Align with a picked reference, then a popup Align from another one.
  const popupAfterDrawer = create({ members: [reference, otherReference, targetA, targetC] });
  const picked = popupAfterDrawer.actions.drawerReferenceOptions('og-1')
    .find(({ anchor: a }) => a.biologicalFeatureId === 'clicked');
  assert.equal(popupAfterDrawer.actions.setDrawerReference('og-1', picked.key), true);
  await popupAfterDrawer.actions.startFromDrawer({ groupId: 'og-1', mode: 'review' });
  await popupAfterDrawer.actions.applyDraft();
  await startReview(popupAfterDrawer, otherReference);
  await popupAfterDrawer.actions.applyDraft();
  assert.equal(popupAfterDrawer.state.similarityAlignmentPlan.value.reference.biologicalFeatureId, 'other');
  assert.equal(popupAfterDrawer.actions.drawerReferenceKey.value, planKey(popupAfterDrawer));

  // Loading a Session installs its own plan.
  const sessionLoad = create({ members: [reference, otherReference, targetA, targetC] });
  await startReview(sessionLoad);
  await sessionLoad.actions.applyDraft();
  sessionLoad.actions.setDrawerReference('og-1', '');
  flushWatchers();
  // The Session saved another drawer group; its plan selects the plan's group.
  sessionLoad.state.selectedOrthogroupId.value = 'og-saved';
  const loaded = JSON.parse(JSON.stringify(sessionLoad.state.similarityAlignmentPlan.value));
  sessionLoad.state.similarityAlignmentPlan.value = loaded;
  flushWatchers();
  assert.equal(sessionLoad.state.selectedOrthogroupId.value, 'og-1');
  assert.equal(sessionLoad.actions.drawerReferenceKey.value, planKey(sessionLoad));
  // A group picked after that plan stays until the next plan.
  sessionLoad.state.selectedOrthogroupId.value = 'og-saved';
  flushWatchers();
  assert.equal(sessionLoad.state.selectedOrthogroupId.value, 'og-saved');

  // The repair route empties the pick; the next committed plan restores the offer.
  const repaired = create();
  await startReview(repaired);
  await repaired.actions.applyDraft();
  const removed = repaired.currentGroup.members.splice(0, 1);
  assert.equal((await repaired.actions.validateBeforeGenerate()).reason, 'stale-reference');
  repaired.actions.setDrawerReference('og-1', '');
  repaired.currentGroup.members.unshift(...removed);
  await startReview(repaired);
  await repaired.actions.applyDraft();
  assert.equal(repaired.actions.drawerReferenceKey.value, planKey(repaired));
  assert.equal(repaired.actions.drawerDisabledReason('og-1'), '');
});

test('malformed Python projection and initial Worker errors leave the prior Result intact', async () => {
  for (const helper of [
    (_operation, { request }) => ({ result: { ...responseFor(request), schema: 1 } }),
    (_operation, { request }) => {
      const response = responseFor(request);
      response.records[1].candidates[0].orientation = {};
      return { result: response };
    },
    (_operation, { request }) => {
      const response = responseFor(request);
      response.records[1].orientationPolicy = 'preserve';
      return { result: response };
    },
    (_operation, { request }) => {
      const response = responseFor(request);
      response.plan.records[1].effectiveReverseComplement = false;
      return { result: response };
    },
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

test('manual Reverse keeps the plan while invalidating edits clear it with a notice', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const plan = fixture.state.similarityAlignmentPlan.value;
  fixture.state.activeDrawing().linearSeqs[1].region_reverse = true;
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.state.activeDrawing().linearSeqs[1].region_reverse, true);
  assert.equal(fixture.actions.notice.value, '');
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.clearForMutation('record crop changed.', fixture.state.activeDrawingId.value), true);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
  assert.match(fixture.actions.notice.value, /record crop changed/);
  // OV-382: an Undo that restores the plan retires the notice.
  const restored = { ...plan };
  fixture.state.similarityAlignmentPlan.value = restored;
  assert.equal(fixture.actions.notice.value, '');
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
  // OV-342: nothing was clicked, so the notice names the missing reference.
  assert.equal(fixture.actions.notice.value,
    'Alignment needs repair: The alignment reference feature is no longer found in its record.');
});

// OV-3676: the notice is worded where it shows, and the repair note keeps its reason as error text parts that the
// drawer words where it shows, so a UI language change rewords both.
test('the repair notice rewords after a UI language change, and the note keeps its English reason', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.currentGroup.members.splice(0, 1);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'blocked', reason: 'stale-reference' });
  const reason = 'The alignment reference feature is no longer found in its record.';
  await setUiLocale('ja', async () => ({ default: {
    'Alignment needs repair: {reason}': '整列の修復が必要です：{reason}', [reason]: '参照フィーチャーが見つかりません。'
  } }));
  try {
    assert.equal(fixture.actions.notice.value, '整列の修復が必要です：参照フィーチャーが見つかりません。');
    assert.deepEqual(fixture.actions.repair.value.reason, [{ text: reason }]);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// OV-3676: the review's candidate line keeps what its words are made of (the English keys of the evidence kinds, the
// strand and the recommendation reason), so the panel words them in the language it shows, also after a change.
test('the review draft keeps the candidate line\'s English keys when it opens in another language', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => {
      const response = responseFor(request, { ambiguous: true });
      response.records[1].candidates[0].directEvidence = ['rbh', 'helper_kind'];
      return { result: response };
    } });
  await setUiLocale('ja', async () => ({ default: {
    'reciprocal best hit': '相互ベストヒット',
    'the only group representative among the candidates': 'グループの代表である唯一の候補',
    unknown: '不明', Unchanged: '変更なし'
  } }));
  try {
    await startReview(fixture);
    const [row] = fixture.actions.draft.value.rows;
    assert.deepEqual(row.candidates[0].evidence, ['reciprocal best hit', 'helper_kind']);
    assert.equal(row.recommendationReasonLabel, 'the only group representative among the candidates');
    assert.equal(Object.hasOwn(row, 'reasonLabel'), false, 'no stored label that nothing shows');
    // The candidate line and the recommendation as index.html words them.
    const facts = templateBinding(/<span class="alignment-review-facts">\{\{ (.+?) \}\}/);
    const strand = templateBinding(/<span class="alignment-review-facts">.*? · \{\{ candidate\.coordinates \}\} \(\{\{ (.+?) \}\}\)<\/span>/);
    const recommended = templateBinding(/<template v-if="record\.recommendedKey">\{\{ (.+?) \}\}<\/template>/);
    const [candidate] = row.candidates;
    assert.equal(facts({ candidate, t, tm }), '相互ベストヒット, helper_kind');
    assert.equal(strand({ candidate, tm }), candidate.displayedStrand);
    assert.match(recommended({ record: row, t, tm, alignmentFeatureName: reviewView.alignmentFeatureName }), /グループの代表である唯一の候補/);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// OV-3676: a record with no name of its own is named by its position ("Record {n}", also beside two records with one
// name), and a feature without a type "Feature". The draft keeps the position and the missing type as data, so the
// panel, the direction preview and the figure overlay word them in the shown language, also after a change.
test('the review words record and feature name fallbacks in the language shown, also after a change', async () => {
  const untyped = { ...targetC, gene: undefined, proteinId: undefined, sourceProteinId: undefined, type: undefined };
  const fixture = create({ members: [reference, targetA, untyped] });
  // Record a has no name of its own; b and c have one definition.
  fixture.state.activeDrawing().linearSeqs[1].definition = 'Genome B';
  fixture.state.activeDrawing().linearSeqs[2].definition = 'Genome B';
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
  const draft = fixture.actions.draft.value;
  // The overlay is drawn while the review opens in English.
  const overlay = reviewView.alignmentOverlayRecords(draft.rows);
  await setUiLocale('ja', async () => ({ default: {
    'Record {n}': 'レコード{n}', Feature: 'フィーチャー',
    'Select {label}, {coordinates} for {record}': '{record}の{label}（{coordinates}）を選択'
  } }));
  try {
    const names = { alignmentRecordName: reviewView.alignmentRecordName, alignmentFeatureName: reviewView.alignmentFeatureName };
    // The card's record name and its accessible name, as index.html words them.
    const legend = templateBinding(/<legend v-if="card\.kind !== 'reference'" class="sr-only">\{\{ (.+?) \}\}<\/legend>/);
    const shownName = templateBinding(/<span dir="auto" class="alignment-review-record">\{\{ (.+?) \}\}<\/span>/);
    const single = templateBinding(/^\s*<b dir="auto">\{\{ (.*?card\.candidates\[0\].*?) \}\}<\/b>$/m);
    assert.deepEqual(draft.rows.map((record) => legend({ record, ...names })), ['Genome B · レコード2', 'Genome B · レコード3']);
    assert.equal(shownName({ record: draft.reference, ...names }), 'レコード1');
    const untypedCandidate = draft.rows[1].candidates[0];
    assert.equal(single({ card: { candidates: [untypedCandidate] }, ...names }), 'フィーチャー 41..60 bp');
    // The direction preview, a computed of the draft.
    assert.deepEqual(fixture.actions.directionPreview.value.records.map(({ label }) => label),
      ['レコード1', 'Genome B · レコード2', 'Genome B · レコード3']);
    // The overlay's ring names, worded where the overlay writes them.
    assert.equal(overlay[1].candidates[0].accessibleName(), 'Genome B · レコード3のフィーチャー 41..60 bp（41..60 bp）を選択');
    assert.equal(overlay[1].candidates[0].name(), 'フィーチャー 41..60 bp');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// OV-3682: a repair whose helper run fails keeps the normalized parts of the failure, never its raw text, and the
// drawer and the notice word them in the shown language.
test('a failed repair run keeps the failure\'s normalized parts in the repair note', async () => {
  let fail = false;
  const fixture = create({ helper: (_operation, { request }) => {
    if (fail) throw new Error('PRIVATE helper text');
    return { result: responseFor(request) };
  } });
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.currentGroup.members.splice(1, 1, targetB);
  fail = true;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'blocked', reason: 'stale-reference' });
  const { reason } = fixture.actions.repair.value;
  assert.ok(Array.isArray(reason), 'the note keeps error text parts');
  assert.equal(JSON.stringify(reason).includes('PRIVATE'), false, 'no raw exception text');
  assert.equal(fixture.actions.notice.value.includes('PRIVATE'), false);
  const english = errorText(reason);
  assert.match(english, /^The exact reference is no longer available\. /);
  await setUiLocale('ja', async () => ({ default: {
    'Alignment needs repair: {reason}': '整列の修復が必要です：{reason}',
    'The exact reference is no longer available.': '参照フィーチャーは利用できなくなりました。'
  } }));
  try {
    assert.ok(errorText(reason).startsWith('参照フィーチャーは利用できなくなりました。'), errorText(reason));
    const note = templateBinding(/<p class="mt-1">\{\{ (.+?similarityAlignmentRepair\.reason.*?) \}\}<\/p>/);
    assert.equal(note({ similarityAlignmentRepair: fixture.actions.repair.value, errorText, tm }), errorText(reason));
    assert.ok(fixture.actions.notice.value.startsWith('整列の修復が必要です：参照フィーチャーは利用できなくなりました。'),
      fixture.actions.notice.value);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// E1: the plan is Linear draft state; a Circular Generate does not draw it, so
// a plan whose reference left the Circular Result's groups does not block it.
test('a Circular Generate is not blocked by the Linear alignment plan', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const plan = fixture.state.similarityAlignmentPlan.value;
  fixture.currentGroup.members.splice(0, 1);
  fixture.state.mode.value = 'circular';
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  fixture.state.mode.value = 'linear';
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
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
  // A target feature that left its record still blocks (OV-341).
  fixture.state.featureCatalog.value.items[0].biologicalFeatures = [reference, targetC];
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.actions.notice.value,
    'Alignment needs repair: An aligned target feature is no longer found in its record.');
});

// OV-341 (Owner decision D-46 = A): a plan names anchor features, not group
// membership. A comparison setting that regroups keeps the plan while every
// anchor feature still resolves in the shown Result.
test('a regrouped Similarity Group keeps the plan while its anchors resolve (OV-341)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const plan = fixture.state.similarityAlignmentPlan.value;
  const calls = fixture.helperCalls.length;
  // The reference left the plan's group; its feature is still in record a.
  fixture.currentGroup.members = [targetA, targetC];
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  // A no-edit Generate after that one succeeds as well.
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  // A target left the group instead.
  fixture.currentGroup.members = [reference, targetC];
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.helperCalls.length, calls);
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  assert.equal(fixture.actions.repair.value, null);
  assert.equal(fixture.actions.notice.value, '');

  // A plan that the shown Result did not draw, whose group left the catalog.
  const loaded = create({ getCommittedRequest: () => renderRequest });
  await startReview(loaded);
  await loaded.actions.applyDraft();
  loaded.state.similarityAlignmentPlan.value = structuredClone(plan);
  loaded.currentGroup.id = 'no-longer-in-catalog';
  assert.deepEqual(await loaded.actions.validateBeforeGenerate(), { status: 'ok' });
});

// An older Session's Result can have no feature catalog, and its committed
// request need not carry the plan. The client cannot check anchors then, so
// the typed render validates them (OV-341, #997 CI).
test('without a feature catalog the plan is not blocked as stale', async () => {
  const fixture = create({ getCommittedRequest: () => renderRequest });
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.state.featureCatalog.value = null;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  fixture.currentGroup.id = 'no-longer-in-catalog';
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' });
  assert.equal(fixture.actions.repair.value, null);
});

test('a gone target next to a regrouped target blocks without a helper call (OV-342)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const calls = fixture.helperCalls.length;
  // target-a left record b; target-c is still in record c but left the group.
  fixture.currentGroup.members = [reference];
  fixture.state.featureCatalog.value.items[0].biologicalFeatures = [reference, targetC];
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.actions.notice.value,
    'Alignment needs repair: An aligned target feature is no longer found in its record.');
  assert.equal(fixture.helperCalls.length, calls);
});

test('a regrouped plan whose reference feature is gone names it in the repair notice (OV-342)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  const previous = fixture.state.results.value;
  fixture.currentGroup.members = [targetA, targetC];
  fixture.state.featureCatalog.value.items[0].biologicalFeatures = [targetA, targetC];
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.actions.notice.value,
    'Alignment needs repair: The alignment reference feature is no longer found in its record.');
  assert.equal(fixture.actions.repair.value.kind, 'reference');
  assert.equal(fixture.state.results.value, previous);
});

test('candidate labels retain biological metadata and displayed record names', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  fixture.state.activeDrawing().linearSeqs[1].definition = 'Genome B';
  await startReview(fixture);
  const row = fixture.actions.draft.value.rows[0];
  assert.equal(reviewView.alignmentRecordName(row.recordNames).label, 'Genome B');
  assert.equal(reviewView.alignmentFeatureName(row.candidates[0].feature).label, 'gene-target-a · target-a');
  assert.equal(row.candidates[0].coordinates, '21..40 bp');
  assert.equal(row.candidates[0].displayedStrand, '-');
});


test('exclusive modes, reference Custom, Select and Skip are local with zero Worker', async () => {
  const f=create({members:[reference,targetA,targetB,targetC]});await startReview(f);
  f.actions.setDirectionMode('left');
  assert.deepEqual(f.actions.directionPreview.value.records.map(r=>r.afterReverseComplement),[true,false,false]);
  f.actions.setDirectionMode('custom');f.actions.setCustomDirection('a','left');f.actions.setCustomDirection('b','right');
  assert.deepEqual(f.actions.directionPreview.value.records.map(r=>r.afterReverseComplement),[true,true,false]);
  f.actions.skipRecord('b');assert.equal(f.actions.directionPreview.value.records[1].afterReverseComplement,false);
  f.actions.setDirectionMode('keep');assert.deepEqual(f.actions.draft.value.intent,{mode:'keep'});
  assert.equal(f.helperCalls.length,1);assert.equal(f.generationCalls.length,0);
  // An unknown-strand anchor still moves its record; only its direction stays (PD-OI-027).
  const unknown=f.actions.directionPreview.value.records.find(r=>r.exclusion==='unknown_strand');
  assert.match(unknown.exclusionLabel,/^Direction unchanged: /);
  // Each count line is its own sentence in the joined effect, never run together.
  assert.match(f.actions.directionPreview.value.effect,/\. 1 record unchanged\. 1 record keeps its direction\.$/);
});
test('Apply and both Reset scopes own one History transaction and consume receipt',async()=>{
  for(const scope of ['positions','positions-and-directions']){
    const f=create();await startReview(f);f.actions.setDirectionMode('left');await f.actions.applyDraft();
    assert.equal(f.history.getUndoCount(),1);assert.equal(f.state.activeDrawing().linearSeqs[0].region_reverse,true);
    const receipt=structuredClone(f.state.similarityAlignmentResetReceipt.value);assert.deepEqual(receipt.directions,[{recordKey:'a',before:false,after:true}]);
    assert.ok(receipt.referenceDeltaX);const calls=f.helperCalls.length;
    await f.actions.resetAlignment(scope);assert.equal(f.helperCalls.length,calls+Number(scope==='positions-and-directions'));
    assert.equal(f.state.similarityAlignmentResetReceipt.value,null);assert.equal(f.state.similarityAlignmentPlan.value,null);
    assert.equal(f.state.activeDrawing().linearSeqs[0].region_reverse,scope==='positions');
    await f.history.undo();assert.deepEqual(f.state.similarityAlignmentResetReceipt.value,receipt);
    assert.equal(f.state.activeDrawing().linearSeqs[0].region_reverse,true);await f.history.redo();assert.equal(f.state.similarityAlignmentResetReceipt.value,null);
  }
});
test('changed final reference correction refreshes preview and requires another Apply',async()=>{
  const f=create({mutateProjection:(response,count)=>{
    if(count>=2)response.projection.records[0].variants[1].anchors[0].centerX+=5;
  }});await startReview(f);f.actions.setDirectionMode('left');
  const first=f.actions.directionPreview.value.reference.deltaX;
  assert.equal((await f.actions.applyDraft()).status,'reviewing');assert.equal(f.generationCalls.length,0);
  assert.equal(f.helperCalls.length,2);assert.equal(f.actions.directionPreview.value.reference.deltaX,first-5);
  assert.equal((await f.actions.applyDraft()).status,'ok');assert.equal(f.helperCalls.length,3);assert.equal(f.generationCalls.length,1);
});
test('a changed final source binding rejects corrections without committing anything',async()=>{
  const f=create({mutateProjection:(response,count)=>{if(count===2)response.projection.binding='b'.repeat(64);}});
  await startReview(f);f.actions.setDirectionMode('right');assert.equal((await f.actions.applyDraft()).status,'stale');
  assert.equal(f.history.getUndoCount(),0);assert.equal(f.generationCalls.length,0);assert.equal(f.state.similarityAlignmentResetReceipt.value,null);
});
test('empty modern evidence and missing historical evidence disable combined Reset with different reasons',async()=>{
  const f=create();await startAlign(f);assert.match(f.actions.resetPreview.value.disabledReason,/did not change any record direction/);
  assert.equal(f.actions.resetPreview.value.directionsLabel,'Positions and directions');
  f.state.similarityAlignmentResetReceipt.value=null;assert.match(f.actions.resetPreview.value.disabledReason,/does not record which directions/);
  assert.equal((await f.actions.resetAlignment('positions-and-directions')).status,'rejected');
  assert.equal((await f.actions.resetAlignment()).status,'ok');
});
test('canceled and stale direction generations leave orientations, receipt and History unchanged', async () => {
  for (const outcome of ['canceled', 'stale']) {
    const pending = deferred();const fixture=create({generation:()=>pending.promise});await startReview(fixture);
    fixture.actions.setDirectionMode('right');const operation=fixture.actions.applyDraft();await new Promise(resolve=>setImmediate(resolve));
    assert.equal(fixture.state.activeDrawing().linearSeqs[1].region_reverse,false);if(outcome==='canceled')fixture.actions.cancel();pending.resolve({status:outcome});
    assert.equal((await operation).status,outcome==='canceled'?'stale':outcome);assert.equal(fixture.history.getUndoCount(),0);
    assert.equal(fixture.state.similarityAlignmentResetReceipt.value,null);
  }
});

test('Align A then B resets to B before, including reference and only affected later edits', async () => {
  const f = create();
  await startReview(f);f.actions.setDirectionMode('left');await f.actions.applyDraft();
  const afterA = f.state.activeDrawing().linearSeqs.map(row => row.region_reverse);
  await startReview(f);f.actions.setDirectionMode('right');
  const firstB=await f.actions.applyDraft();
  if(firstB.status==='reviewing'){assert.equal(f.history.getUndoCount(),1);assert.equal((await f.actions.applyDraft()).status,'ok');}
  else assert.equal(firstB.status,'ok');
  const receiptB = structuredClone(f.state.similarityAlignmentResetReceipt.value);
  assert.deepEqual(receiptB.directions, [
    {recordKey:'a',before:true,after:false}, {recordKey:'b',before:false,after:true}
  ]);
  assert.equal(f.history.getUndoCount(),2);
  // A later affected edit is replaced, while a skipped/unknown record's edit survives.
  f.state.activeDrawing().linearSeqs[0].region_reverse=true;
  f.state.activeDrawing().linearSeqs[2].region_reverse=true;
  // Model the ordinary Generate owner committing those manual directions.
  f.generationCalls.at(-1).canonical.renderRequest.records[0].presentation.reverseComplement=true;
  f.generationCalls.at(-1).canonical.renderRequest.records[2].presentation.reverseComplement=true;
  f.generationCalls.at(-1).canonical.renderRequest.records[0].presentation.label='<i>Exact displayed record</i>';
  f.state.activeDrawing().linearSeqs[0].definition='Source-only name';
  assert.equal(f.actions.resetPreview.value.targets[0].label,'Exact displayed record');
  assert.equal(f.actions.resetPreview.value.targets[0].name,'Exact displayed record');
  assert.equal(f.actions.resetPreview.value.directionsLabel,'Positions and directions (2 records)');
  assert.equal(f.actions.resetPreview.value.targets[0].laterManualEdit,true);
  // The direction words go through the catalogs (the template never holds them).
  const {current,restored,currentWord,restoredWord}=f.actions.resetPreview.value.targets[0];
  assert.deepEqual([currentWord,restoredWord],[current,restored].map(reverse=>reverse?'reverse':'forward'));
  assert.equal((await f.actions.resetAlignment('positions-and-directions')).status,'ok');
  assert.deepEqual(f.state.activeDrawing().linearSeqs.map(row=>row.region_reverse),[...afterA.slice(0,2),true]);
  assert.equal(f.state.similarityAlignmentResetReceipt.value,null);
  assert.equal(f.history.getUndoCount(),3);
  await f.history.undo();assert.deepEqual(f.state.similarityAlignmentResetReceipt.value,receiptB);
  await f.history.undo();assert.deepEqual(f.state.activeDrawing().linearSeqs.map(row=>row.region_reverse),afterA);
});

test('failed Reset retains the original receipt and Result for an explicit retry', async () => {
  let fail=false;
  const f=create({generation:async()=>fail?{status:'error',error:{summary:'Reset render rejected.'}}:{status:'ok'}});
  await startReview(f);f.actions.setDirectionMode('left');await f.actions.applyDraft();
  const before=structuredClone({receipt:f.state.similarityAlignmentResetReceipt.value,
    plan:f.state.similarityAlignmentPlan.value,translations:f.state.linearRecordTranslations.value,
    sequences:f.state.activeDrawing().linearSeqs,results:f.state.results.value});
  fail=true;assert.equal((await f.actions.resetAlignment('positions-and-directions')).status,'error');
  assert.deepEqual({receipt:f.state.similarityAlignmentResetReceipt.value,
    plan:f.state.similarityAlignmentPlan.value,translations:f.state.linearRecordTranslations.value,
    sequences:f.state.activeDrawing().linearSeqs,results:f.state.results.value},before);
  assert.equal(f.history.getUndoCount(),1);
  fail=false;assert.equal((await f.actions.resetAlignment('positions-and-directions')).status,'ok');
  assert.equal(f.history.getUndoCount(),2);assert.equal(f.state.similarityAlignmentResetReceipt.value,null);
});

test('combined Reset uses current geometry for y without storing a translation snapshot', async () => {
  const f=create({mutateProjection:(response,count)=>{
    if(count===3)response.projection.records[0].variants[0].axisY+=20;
  }});
  await startReview(f);f.actions.setDirectionMode('left');await f.actions.applyDraft();
  const prior=f.state.linearRecordTranslations.value.find(row=>row.recordKey==='a').y;
  assert.equal((await f.actions.resetAlignment('positions-and-directions')).status,'ok');
  assert.equal(f.state.linearRecordTranslations.value.find(row=>row.recordKey==='a').y,prior-20);
  assert.equal(f.helperCalls.length,3);
});

test('stable reorder retains receipt, invalidation consumes it through the existing plan owner', async () => {
  const f=create();await startReview(f);f.actions.setDirectionMode('left');await f.actions.applyDraft();
  const receipt=structuredClone(f.state.similarityAlignmentResetReceipt.value);
  f.actions.retainForStableReorder(['c','b','a']);
  assert.deepEqual(f.state.similarityAlignmentResetReceipt.value,receipt);
  f.actions.clearForMutation('Source changed.', f.state.activeDrawingId.value);
  assert.equal(f.state.similarityAlignmentPlan.value,null);
  assert.equal(f.state.similarityAlignmentResetReceipt.value,null);
});


test('missing caller diagnostic cannot borrow an unrelated global error',async()=>{
  const fixture=create({generation:async()=>({status:'error'})});
  fixture.state.errorLog.value={code:'INPUT_UNREADABLE',operation:'export-svg',summary:'unrelated'};
  const outcome=await startAlign(fixture);
  assert.equal(outcome.status,'error');
  assert.equal(outcome.error.code,'UNKNOWN');
  assert.equal(outcome.error.operation,'align');
  assert.match(outcome.error.summary,/without recognized diagnostic information/);
  assert.doesNotMatch(outcome.error.summary,/input could not be read|network/i);
  assert.equal(fixture.actions.dialogOpen.value,true);
});

test('Review and the active-plan inspector use names, never internal keys or hashes', async () => {
  const unnamed = { ...targetC, gene: undefined, proteinId: undefined, sourceProteinId: undefined, type: 'CDS' };
  const fixture = create({ members: [reference, targetA, unnamed] });
  Object.assign(fixture.state.activeDrawing().linearSeqs[0], { definition: 'Genome A', accession: 'ACC-A' });
  fixture.state.activeDrawing().linearSeqs[1].definition = 'Genome B';
  fixture.state.activeDrawing().linearSeqs[2].definition = 'Genome B';
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
  const draft = fixture.actions.draft.value;
  assert.equal(reviewView.alignmentRecordName(draft.reference.recordNames).label, 'Genome A · ACC-A');
  assert.deepEqual(draft.rows.map(({ recordNames }) => reviewView.alignmentRecordName(recordNames).label),
    ['Genome B · Record 2', 'Genome B · Record 3']);
  assert.equal(reviewView.alignmentFeatureName(draft.rows[1].candidates[0].feature).label, 'CDS 41..60 bp');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  const inspector = fixture.actions.activePlanInspector.value;
  assert.equal(inspector.reference.recordLabel, 'Genome A · ACC-A');
  assert.equal(inspector.reference.label, 'gene-clicked · clicked');
  assert.deepEqual(inspector.records.map(({ recordLabel, anchorLabel }) => [recordLabel, anchorLabel]), [
    ['Genome A · ACC-A', 'gene-clicked · clicked'],
    ['Genome B · Record 2', 'gene-target-a · target-a'],
    ['Genome B · Record 3', 'CDS 41..60 bp']
  ]);
  // <record-name> lays the parts out: the name shrinks, the detail stays whole.
  assert.deepEqual(inspector.records.map(({ recordName, recordDetail }) => [recordName, recordDetail]), [
    ['Genome A', 'ACC-A'], ['Genome B', 'Record 2'], ['Genome B', 'Record 3']
  ]);
});

// OV-2883: the inspector reads a record's orientation where its File card holds it: on the windows the card lists.
test('the active-plan inspector reads a card drawn in windows as reversed by its windows', async () => {
  const fixture = create({ members: [reference, targetA, targetC] });
  await startReview(fixture);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  const card = fixture.state.activeDrawing().linearSeqs[1];
  const reversed = () => fixture.actions.activePlanInspector.value.records
    .find(({ recordKey }) => recordKey === card.uid).reversedFromSource;
  card.region_reverse = false;
  card.windows = [{ start: 1, end: 10, reverseComplement: true, column: null }];
  assert.equal(reversed(), true);
  card.windows = [{ start: 1, end: 10, reverseComplement: false, column: null }];
  card.region_reverse = true;
  assert.equal(reversed(), false, 'a card that lists windows is not oriented by its Reverse complement');
});

// UX-B06: the drawer reference options name the record and the gene, not the
// record key and the feature hash; the active plan is one line.
test('drawer reference options and the active-plan line use record and gene names', async () => {
  const inparalog = { ...member('a', 'twin', 9, '-'), gene: 'gene-clicked' };
  const fixture = create({ members: [reference, inparalog, targetA, targetC] });
  Object.assign(fixture.state.activeDrawing().linearSeqs[0], { definition: 'Genome A', accession: 'ACC-A' });
  fixture.state.activeDrawing().linearSeqs[1].definition = 'Genome B';
  assert.deepEqual(fixture.actions.drawerReferenceOptions('og-1').map(({ label }) => label), [
    'ACC-A · gene-clicked · 51..70 (+)',
    'ACC-A · gene-clicked · 101..120 (-)',
    'Genome B · gene-target-a',
    'Record 3 · gene-target-c'
  ]);
  await startReview(fixture);
  await fixture.actions.applyDraft();
  assert.equal(fixture.actions.activePlanInspector.value.summary, 'Aligned to gene-clicked in ACC-A');
});

// Last in this file: without the own-key read, `__proto__` handed `Object.prototype` to the draft's deep freeze.
test('an evidence kind outside the label table shows as written, also an Object.prototype name (OV-573)', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => {
      const response = responseFor(request, { ambiguous: true });
      response.records[1].candidates[0].directEvidence = ['constructor', 'rbh', 'related_homolog', '__proto__', ''];
      return { result: response };
    } });
  await startReview(fixture);
  assert.deepEqual(fixture.actions.draft.value.rows[0].candidates[0].evidence,
    ['constructor', 'reciprocal best hit', 'related_homolog', '__proto__', '']);
  assert.equal(Object.isFrozen(Object.prototype), false);
});

// OV-572: a record key is data (a Session's uid), so the Custom select reads the direction
// table by own key, through the alignment owner. The expression is the one index.html binds.
test('the Custom direction select reads its record by own key (OV-572)', async () => {
  const f = create({ members: [reference, targetA, targetB, targetC] });
  await startReview(f);
  f.actions.setDirectionMode('custom');
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const binding = html.match(/:value="([^"]+)"\s+:aria-label="t\('Direction for \{record\}', \{ record: card\.direction\.label \}\)"/)?.[1];
  assert.ok(binding, 'the Custom direction select is found');
  const selectValue = new Function(
    'similarityAlignmentDraft', 'similarityAlignmentCustomDirection', 'card', `return (${binding});`
  );
  // What the composition root hands the template (app/app-setup.js).
  const value = (recordKey) => selectValue(f.actions.draft.value, f.actions.customDirection, { recordKey });
  assert.equal(value('b'), 'keep', 'no choice yet: Keep');
  assert.equal(value('constructor'), 'keep', 'an inherited member names no choice');
  assert.equal(value('__proto__'), 'keep');
  assert.equal(value('toString'), 'keep');
  f.actions.setCustomDirection('b', 'left');
  assert.equal(value('b'), 'left');
  assert.equal(value('a'), 'keep');
  // Outside Custom the draft has no table; the reader still answers.
  f.actions.setDirectionMode('keep');
  assert.equal(f.actions.customDirection('b'), 'keep');
  // The template reads the table through the owner, and the root binds the reader.
  assert.doesNotMatch(html, /byRecordKey\s*\[/);
  const appSetup = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
  assert.match(appSetup, /similarityAlignmentCustomDirection: similarityAlignmentActions\.customDirection,/);
});

// D-38 (PD-OI-039, PD-OI-090): a review belongs to the drawing it was opened in, by its drawing ID.
// Another drawing here shows the same committed facts, as a Duplicate does: the review stays open
// over it, draws nothing on it and applies nothing to it. D-39 (OV-3585): there Apply is disabled
// with a reason that names the drawing to show, and keeps the review and its choices. Shown again,
// the review's own drawing takes the Apply.
const showDrawing = (fixture, drawingId) => {
  fixture.state.activeDrawingId.value = drawingId;
  fixture.actions.arriveAtDrawing();
};
test('a review stays with its own drawing across a switch and applies only there (D-38)', async () => {
  const fixture = create();
  await startReview(fixture);
  fixture.actions.summary.value = { text: 'departing summary' };
  showDrawing(fixture, 'linear-2');
  assert.equal(fixture.actions.status.value, 'reviewing', 'the switch keeps the review open');
  assert.notEqual(fixture.actions.draft.value, null, 'and its choices');
  assert.equal(fixture.actions.summary.value, null, 'the summary was the departing drawing\'s');
  assert.equal(fixture.actions.reviewDrawingId(), 'linear-1');
  assert.equal(fixture.actions.isDraftArtifactCurrent(), false, 'the draft does not show on another drawing');
  showDrawing(fixture, 'linear-1');
  assert.equal(fixture.actions.isDraftArtifactCurrent(), true, 'it shows again on its own drawing');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.generationCalls.length, 1);
  assert.equal(fixture.state.similarityAlignmentPlan.value.groupId, 'og-1');

  await startReview(fixture);
  const plan = fixture.state.similarityAlignmentPlan.value;
  showDrawing(fixture, 'linear-2');
  const choices = fixture.actions.draft.value;
  assert.equal(fixture.actions.canApply.value, false, 'Apply is disabled over another drawing');
  assert.equal(fixture.actions.applyDisabledReason.value, 'Show Drawing linear-1 to apply this alignment.');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'rejected' }, 'Apply over another drawing changes none');
  assert.equal(fixture.generationCalls.length, 1);
  assert.equal(fixture.state.similarityAlignmentPlan.value, plan);
  assert.equal(fixture.actions.draft.value, choices, 'the review keeps its choices');
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.error.value, null, 'and says nothing changed');
  showDrawing(fixture, 'linear-1');
  assert.equal(fixture.actions.canApply.value, true);
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'ok' });
  assert.equal(fixture.generationCalls.length, 2);
});

test('an edit ends a review only when it edits the review\'s drawing (D-38)', async () => {
  const writes = [];
  const fixture = create({ writeArtifactValues: (drawingId, values) => writes.push([drawingId, values]) });
  await startReview(fixture);
  showDrawing(fixture, 'linear-2');
  fixture.actions.clearForMutation('record crop changed.', 'linear-2');
  assert.equal(fixture.actions.status.value, 'reviewing', 'an edit of the shown drawing leaves another drawing\'s review');
  fixture.actions.clearForMutation('record crop changed.', 'linear-3');
  assert.equal(fixture.actions.status.value, 'reviewing', 'and so does an edit of a third drawing');
  fixture.actions.clearForMutation('record crop changed.', 'linear-1');
  assert.equal(fixture.actions.status.value, 'idle', 'an edit of the review\'s drawing ends it');
  assert.equal(fixture.actions.draft.value, null);
  assert.deepEqual(writes.map(([drawingId]) => drawingId), ['linear-3', 'linear-1'], 'each hidden drawing\'s plan clears in its slot');

  // Cancel over another drawing ends the review and changes no drawing.
  showDrawing(fixture, 'linear-1');
  await startReview(fixture);
  showDrawing(fixture, 'linear-2');
  assert.deepEqual(fixture.actions.cancel(), { status: 'canceled' });
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.generationCalls.length, 0);
  assert.equal(fixture.state.similarityAlignmentPlan.value, null);
});

// D-40 (S1): work on the shown drawing never ends, replaces or rewrites another drawing's open review. A Generate of
// another drawing whose plan needs repair draws nothing and says why with the repair note; the open review stays.
test('a Generate of another drawing whose plan needs repair keeps the open review (D-40)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  await startReview(fixture);
  const review = fixture.actions.draft.value;
  showDrawing(fixture, 'linear-2');
  fixture.currentGroup.members.splice(1, 1, targetB);
  const helpers = fixture.helperCalls.length;
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'blocked', reason: 'stale-target' });
  assert.equal(fixture.actions.draft.value, review, 'the other drawing\'s review keeps its draft');
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.reviewDrawingId(), 'linear-1');
  assert.equal(fixture.helperCalls.length, helpers, 'no repair review is resolved for the shown drawing');
  assert.equal(fixture.actions.repair.value.kind, 'targets');
  // D-41 (S3): the notice says why nothing opened, and names the drawing whose review to finish.
  assert.equal(fixture.actions.notice.value,
    'Finish or cancel the alignment review in Drawing linear-1, then Generate this drawing again.');

  fixture.currentGroup.members.splice(0, 1);
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'blocked', reason: 'stale-reference' });
  assert.equal(fixture.actions.draft.value, review, 'a stale reference of the shown drawing keeps it too');
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.repair.value.kind, 'reference');
});

const GENERATING = { status: 'busy', reason: 'Generating diagram. Retry after generation finishes.' };

// Carried (c): a run refuses Align and Review before the resolve, so an ambiguous group opens no review either; the
// drawer's controls are disabled with the run's reason (carried (b)).
test('a run refuses Align and Review before the resolve and the drawer states its reason', async () => {
  const fixture = create({ members: [reference, targetA, targetB, targetC],
    helper: (_operation, { request }) => ({ result: responseFor(request, { ambiguous: true }) }) });
  fixture.state.runBlock.value = GENERATING;
  assert.equal(fixture.actions.drawer.value.disabledReason, GENERATING.reason);
  assert.deepEqual(await startReview(fixture), GENERATING);
  assert.deepEqual(await startAlign(fixture), GENERATING);
  assert.equal(fixture.helperCalls.length, 0, 'no resolve runs');
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.actions.status.value, 'idle');
  assert.equal(fixture.state.errorLog.value, null, 'a refusal is no failure');
  fixture.state.runBlock.value = null;
  assert.deepEqual(await startReview(fixture), { status: 'reviewing' });
});

// OV-2447: a review opened before a Generate stays open; while the run holds it and after the Generate commits
// another Result, its Apply is disabled with the reason, and an Apply that still runs fails with a recognized code.
test('a review across a Generate disables Apply with the reason and a stale Apply fails as ALIGNMENT_STALE', async () => {
  const fixture = create();
  await startReview(fixture);
  assert.equal(fixture.actions.canApply.value, true);
  fixture.state.runBlock.value = GENERATING;
  assert.equal(fixture.actions.canApply.value, false);
  assert.equal(fixture.actions.applyDisabledReason.value, GENERATING.reason);
  assert.deepEqual(await fixture.actions.applyDraft(), GENERATING);
  assert.equal(fixture.helperCalls.length, 1, 'a refused Apply runs no resolve');
  fixture.replaceCommittedResources();
  fixture.state.runBlock.value = null;
  assert.equal(fixture.actions.dialogOpen.value, true, 'the review stays open');
  assert.equal(fixture.actions.canApply.value, false);
  assert.equal(fixture.actions.applyDisabledReason.value,
    'The diagram changed after this review opened. Start alignment again.');
  assert.deepEqual(await fixture.actions.applyDraft(), { status: 'stale' });
  assert.equal(fixture.state.errorLog.value.code, 'ALIGNMENT_STALE');
  assert.equal(errorText(fixture.state.errorLog.value.summaryParts),
    'The diagram changed after this review opened. Start alignment again.');
});

// OV-3849: an Align or Review the user can fix fails with its own code and sentence, not as UNKNOWN.
test('an alignment without one exact reference, of a group not drawn, or not validated fails with its own code', async () => {
  const cases = [
    ['ALIGNMENT_REFERENCE', create(), (fixture) => startAlign(fixture, null)],
    ['ALIGNMENT_GROUP_MISSING', create(), (fixture) => fixture.actions.startFromPopup({ groupId: 'og-gone', reference })],
    ['ALIGNMENT_NOT_VALIDATED', create({ mutateProjection: (response, count) => {
      if (count !== 1) return;
      const orientations = response.projection.geometryOrientations;
      orientations.a = !orientations.a;
    } }), (fixture) => startAlign(fixture)]
  ];
  const outcomes = [];
  /** @type {string[]} */
  const texts = [];
  for (const [, fixture, run] of cases) {
    const { status } = await run(fixture);
    const { code, operation, summaryParts } = fixture.state.errorLog.value;
    outcomes.push({ status, code, operation, drawn: fixture.generationCalls.length, after: fixture.actions.status.value });
    texts.push(errorText(summaryParts));
  }
  // A failed automatic Apply ends the work as the failed resolve does (`start` awaits it).
  assert.deepEqual(outcomes, cases.map(([code]) => ({ status: 'rejected', code, operation: 'align', drawn: 0, after: 'idle' })));
  // Each names the button the user selects instead, as the Similarity groups tab labels it.
  assert.deepEqual(texts.slice(0, 2), [
    'The clicked feature matches no member of this Similarity group, or more than one. Select Align… in the Similarity groups tab after choosing the exact reference record and feature.',
    'This Similarity group is not in the shown diagram. Select Align… in the Similarity groups tab after choosing a group and its exact reference record and feature.'
  ]);
});

// D-40 (S2, N1): a switch waits for the alignment work in flight, with the reason the Apply line gives.
test('the alignment owner names the work a drawing switch waits for (D-40)', async () => {
  let release;
  const fixture = create({ helper: async (_operation, { request }) => {
    await new Promise((resolve) => { release = resolve; });
    return { result: responseFor(request) };
  } });
  assert.equal(fixture.actions.workBusy(), null);
  assert.equal(fixture.actions.busy.value, false);
  const started = startReview(fixture);
  assert.deepEqual(fixture.actions.workBusy(), { status: 'busy', reason: 'Resolving the current alignment.' });
  assert.equal(fixture.actions.busy.value, true, 'busy reads the same work in flight');
  assert.equal(fixture.actions.workReason.value, 'Resolving the current alignment.');
  release();
  await started;
  assert.equal(fixture.actions.workBusy(), null, 'an open review does not hold a switch');
  const applied = fixture.actions.applyDraft();
  assert.deepEqual(fixture.actions.workBusy(), { status: 'busy', reason: 'Applying the alignment plan.' });
  release();
  await applied;
  assert.equal(fixture.actions.workBusy(), null);
});

// D-40 (S3): a removed drawing takes its review with it, as Cancel ends it; another drawing's removal leaves it.
test('removing the review\'s drawing ends its review, another removal leaves it (D-40)', async () => {
  const fixture = create();
  await startReview(fixture);
  showDrawing(fixture, 'linear-2');
  fixture.actions.endReviewOf('linear-3');
  assert.equal(fixture.actions.status.value, 'reviewing');
  fixture.actions.endReviewOf('linear-1');
  assert.equal(fixture.actions.status.value, 'idle');
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.generationCalls.length, 0);
});

// D-40 (N2): a repair review's note is its drawing's. A switch hides it; back on that drawing it shows again, and a
// Cancel there leaves it, as without a switch.
test('a repair review\'s note shows again with its drawing and stays after Cancel (D-40)', async () => {
  const fixture = create();
  await startReview(fixture);
  await fixture.actions.applyDraft();
  fixture.currentGroup.members.splice(1, 1, targetB);
  assert.equal((await fixture.actions.validateBeforeGenerate()).reason, 'stale-target');
  assert.equal(fixture.actions.repair.value.kind, 'targets');
  showDrawing(fixture, 'linear-2');
  assert.equal(fixture.actions.repair.value, null, 'the note was the departing drawing\'s');
  assert.equal(fixture.actions.notice.value, '');
  showDrawing(fixture, 'linear-1');
  assert.equal(fixture.actions.repair.value.kind, 'targets');
  assert.match(fixture.actions.notice.value, /^Alignment needs repair:/);
  fixture.actions.cancel();
  assert.equal(fixture.actions.repair.value.kind, 'targets', 'Cancel leaves the note');
});

// D-41 (B2, OV-3619): a Generate while a review resolves waits for it, and the review opens when the resolve returns;
// before, the Generate's check made the resolve stale and left 'resolving' for good, which held every drawing switch.
test('a Generate while a review resolves waits, and the review still opens (OV-3619)', async () => {
  let release;
  const fixture = create({ helper: async (_operation, { request }) => {
    await new Promise((resolve) => { release = resolve; });
    return { result: responseFor(request) };
  } });
  const started = startReview(fixture);
  assert.equal(fixture.actions.status.value, 'resolving');
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(),
    { status: 'busy', reason: 'Resolving the current alignment.' }, 'Generate waits with the existing reason');
  release();
  assert.deepEqual(await started, { status: 'reviewing' }, 'the resolve is not made stale');
  assert.equal(fixture.actions.status.value, 'reviewing');
  assert.equal(fixture.actions.dialogOpen.value, true);
  assert.equal(fixture.actions.workBusy(), null, 'nothing holds a switch any more');
  assert.deepEqual(await fixture.actions.validateBeforeGenerate(), { status: 'ok' }, 'an open review does not hold Generate');
  assert.equal(fixture.actions.status.value, 'reviewing', 'and the check leaves the review as it is');
});

// D-41 (B1, OV-3618): a Load's commit ends the review of a drawing it did not keep, as Cancel ends it; a Load that
// keeps the review's drawing ID keeps the review (PD-OI-039 r5). A drawing without a name leaves the Apply line
// without the clause that names it.
test('a Load that drops the review\'s drawing ends the review, one that keeps it keeps it (OV-3618)', async () => {
  const fixture = create({ drawingNameOf: (drawingId) => (drawingId === 'linear-1' ? null : `Drawing ${drawingId}`) });
  await startReview(fixture);
  showDrawing(fixture, 'linear-2');
  assert.doesNotThrow(() => fixture.actions.applyDisabledReason.value);
  assert.doesNotMatch(fixture.actions.applyDisabledReason.value, /^Show /, 'a drawing without a tab is not named');
  fixture.actions.endReviewOutside(['linear-1', 'linear-2']);
  assert.equal(fixture.actions.status.value, 'reviewing', 'the review\'s drawing is kept');
  fixture.actions.endReviewOutside(['linear-2']);
  assert.equal(fixture.actions.status.value, 'idle', 'the review leaves with its drawing');
  assert.equal(fixture.actions.draft.value, null);
  assert.equal(fixture.actions.dialogOpen.value, false);
  assert.equal(fixture.generationCalls.length, 0);
  fixture.actions.endReviewOutside([]);
  assert.equal(fixture.actions.status.value, 'idle', 'without a review nothing ends');
});
