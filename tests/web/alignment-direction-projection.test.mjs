import assert from 'node:assert/strict';
import test from 'node:test';
globalThis.window = { Vue: { ref: () => {}, computed: () => {} } };
const { validateSimilarityAlignmentResolution, projectSimilarityAlignmentDirections } = await import(
  '../../gbdraw/web/js/app/similarity-alignment.js'
);
const anchor = (recordKey) => ({ recordKey, biologicalFeatureId: `${recordKey}-protein`,
  sourceFeatureIndex: 0, stableFeatureSvgId: `stable-${recordKey}` });
const fixture = () => {
  const request = { groupId: 'group', reference: anchor('ref'),
    records: ['ref', 'target'].map((recordKey) => ({ recordKey, region: null,
      presentation: { reverseComplement: false } })),
    members: ['ref', 'target'].map((key) => ({ anchor: anchor(key) })) };
  const records = ['ref', 'target'].map((recordKey, index) => ({ kind: 'decision', recordKey,
    status: index ? 'aligned' : 'reference', rationale: index ? 'only_usable_candidate' : 'reference',
    reviewReason: index ? 'only_usable_candidate' : null, anchor: anchor(recordKey), candidates: index ? [{
      anchor: anchor(recordKey), displayName: 'Target', sourceStart: 0, sourceEnd: 40,
      displayCenter: 20, displayedStrand: 1, hidden: false, representative: false,
      role: '', usable: true, directEvidence: [], strandRelation: 'opposite'
    }] : [] }));
  const response = { schema: 2, status: 'resolved', groupId: 'group', reference: anchor('ref'),
    referenceDisplayedStrand: -1, referenceDisplayCenter: 80, records,
    plan: { schema: 2, groupId: 'group', reference: anchor('ref'), records: records.map(
      ({ recordKey, status, rationale, anchor }) => ({ recordKey, status, rationale, anchor })
    ) }, projection: { binding: 'a'.repeat(64), geometryOrientations: { ref: true, target: true },
      records: ['ref', 'target'].map((recordKey, index) => ({ recordKey, beforeReverseComplement: false,
        base: { x: 30 + index * 20, y: 7 + index * 3 }, beforeAxisY: 100 + index * 100,
        beforeAnchors: [{ anchor: anchor(recordKey), centerX: index ? 40 : 160 }],
        variants: [false, true].map((reverseComplement) => ({ reverseComplement,
          axisY: 100 + index * 100 + (reverseComplement ? 12 : 0), anchors: [{ anchor: anchor(recordKey),
            displayCenter: reverseComplement ? (index ? 80 : 20) : (index ? 20 : 80),
            centerX: reverseComplement ? (index ? 160 : 40) : (index ? 40 : 160),
            displayedStrand: (index ? 1 : -1) * (reverseComplement ? -1 : 1) }] })) })) } };
  return { request, response };
};
const project = (response, request, intent, expectedBinding = response.projection.binding) => (
  projectSimilarityAlignmentDirections({ resolution: validateSimilarityAlignmentResolution(response, request), intent, expectedBinding })
);

test('minority reference alone flips right; reference center and every logical y stay fixed', () => {
  const { response, request } = fixture();
  const result = project(response, request, { mode: 'right' });
  assert.deepEqual(result.records.map(({ afterArrow }) => afterArrow), [1, 1]);
  assert.deepEqual(result.records.map(({ afterReverseComplement }) => afterReverseComplement), [true, false]);
  assert.deepEqual(result.reference, { anchor: anchor('ref'), beforeX: 190, afterX: 190, deltaX: 120 });
  assert.deepEqual(result.records.map(({ translation }) => translation.y), [-5, 10]);
  assert.ok(Object.isFrozen(result.records[0].translation));
});

test('left flips only the positive source target, with no reference correction', () => {
  const { response, request } = fixture();
  const result = project(response, request, { mode: 'left' });
  assert.deepEqual(result.records.map(({ afterReverseComplement }) => afterReverseComplement), [false, true]);
  assert.equal(result.reference.deltaX, 0);
});

test('Custom is exclusive, missing rows Keep, and unknown rows cannot be selected', () => {
  const { response, request } = fixture();
  const result = project(response, request, { mode: 'custom', byRecordKey: { ref: 'right' } });
  assert.deepEqual(result.records.map(({ afterReverseComplement }) => afterReverseComplement), [true, false]);
  assert.throws(() => project(response, request, { mode: 'right', byRecordKey: { target: 'left' } }), /invalid fields/);
  assert.throws(() => project(response, request, { mode: 'custom', byRecordKey: { absent: 'right' } }), /Custom/);
  assert.throws(() => project(response, request, { mode: 'custom', byRecordKey: { target: 'reverse' } }), /Custom/);
});

test('binding changes are explicit stale outcomes and cannot produce corrections', () => {
  const { response, request } = fixture();
  const result = project(response, request, { mode: 'right' }, 'b'.repeat(64));
  assert.deepEqual(result, { status: 'stale', reason: 'binding_changed', binding: 'a'.repeat(64) });
  assert.throws(() => projectSimilarityAlignmentDirections({ resolution: response, intent: { mode: 'keep' } }), /explicit/);
});

test('same material output has the same signature across mode and display-name edits', () => {
  const { response, request } = fixture();
  const right = project(response, request, { mode: 'right' });
  response.records[1].candidates[0].displayName = 'A much longer cosmetic name';
  const custom = project(response, request, { mode: 'custom', byRecordKey: { ref: 'right', target: 'keep' } });
  assert.equal(custom.signature, right.signature);
  assert.notEqual(project(response, request, { mode: 'keep' }).signature, right.signature);
});

test('final changed placement changes signature; validation requires exact final orientation coverage', () => {
  const { response, request } = fixture();
  const preview = project(response, request, { mode: 'right' });
  assert.equal(preview.geometryValidated, false);
  const final = structuredClone(response);
  final.projection.geometryOrientations = { ref: true, target: false };
  final.projection.records[0].variants[1].anchors[0].centerX = 45;
  const validated = project(final, request, { mode: 'right' });
  assert.equal(validated.geometryValidated, true);
  assert.notEqual(validated.signature, preview.signature);
  assert.equal(validated.reference.afterX, preview.reference.beforeX);
});

for (const [description, mutate] of [
  ['inconsistent current source strand', (r) => { r.projection.records[1].variants[0].anchors[0].displayedStrand = -1; }],
  ['nonfinite base', (r) => { r.projection.records[0].base.x = Infinity; }],
  ['unknown record', (r) => { r.projection.records[0].recordKey = 'unknown'; }],
  ['duplicate variant', (r) => { r.projection.records[0].variants[1].reverseComplement = false; }],
  ['duplicate anchor', (r) => { r.projection.records[0].variants[0].anchors.push(r.projection.records[0].variants[0].anchors[0]); }],
  ['missing source anchor', (r) => { r.projection.records[0].variants[1].anchors = []; }],
  ['unknown facts field', (r) => { r.projection.policy = 'right'; }],
  ['inconsistent null center', (r) => { r.projection.records[0].variants[0].anchors[0].displayCenter = null; }]
]) test(`strict shared admission rejects ${description}`, () => {
  const { response, request } = fixture();
  mutate(response);
  assert.throws(() => validateSimilarityAlignmentResolution(response, request));
});


test('local Select and Skip share cached Python eligibility and need no Worker', () => {
  const { request, response } = fixture();
  const alternative = { ...anchor('target'), biologicalFeatureId: 'other', sourceFeatureIndex: 1, stableFeatureSvgId: 'other-stable' };
  request.members.push({ anchor: alternative });
  const candidate = { ...response.records[1].candidates[0], anchor: alternative, displayedStrand: -1,
    displayCenter: 30, strandRelation: 'same', displayName: 'Alternative' };
  response.records[1] = { kind: 'ambiguous', recordKey: 'target',
    candidates: [response.records[1].candidates[0], candidate], directRbhCandidates: [],
    recommendedAnchor: anchor('target'), recommendationReason: 'deterministic_candidate_1' };
  response.status = 'ambiguous'; response.plan = null;
  const facts = response.projection.records[1];
  facts.beforeAnchors.push({ anchor: alternative, centerX: 60 });
  facts.variants[0].anchors.push({ anchor: alternative, displayCenter: 30, centerX: 60, displayedStrand: -1 });
  facts.variants[1].anchors.push({ anchor: alternative, displayCenter: 70, centerX: 140, displayedStrand: 1 });
  const resolution = validateSimilarityAlignmentResolution(response, request);
  const run = (choices) => projectSimilarityAlignmentDirections({ resolution, expectedBinding: resolution.projection.binding,
    intent: { mode: 'right' }, choices });
  const unresolved = run([]);
  assert.equal(unresolved.selectionComplete, false);
  assert.equal(unresolved.records[1].exclusion, 'selection_required');
  const selected = run([{ recordKey: 'target', kind: 'select', anchor: alternative }]);
  assert.equal(selected.selectionComplete, true);
  assert.equal(selected.records[1].afterReverseComplement, true);
  assert.equal(selected.records[1].afterArrow, 1);
  assert.deepEqual(selected.records[1].anchor, alternative);
  const skipped = run([{ recordKey: 'target', kind: 'skip', anchor: null }]);
  assert.equal(skipped.records[1].afterReverseComplement, false);
  assert.equal(skipped.records[1].exclusion, 'skipped_by_user');
  assert.notEqual(selected.signature, skipped.signature);
  assert.throws(() => run([{ recordKey: 'ref', kind: 'skip', anchor: null }]), /Invalid local/);
  assert.throws(() => run([{ recordKey: 'target', kind: 'select', anchor: { ...alternative, sourceFeatureIndex: 8 } }]), /usable Python/);
});


test('a changed logical reference destination requires refreshed preview even when delta is unchanged', () => {
  const { response, request } = fixture();
  const preview = project(response, request, { mode: 'right' });
  const refreshed = structuredClone(response);
  refreshed.projection.records[0].beforeAnchors[0].centerX += 5;
  refreshed.projection.records[0].variants.forEach((variant) => { variant.anchors[0].centerX += 5; });
  const final = project(refreshed, request, { mode: 'right' });
  assert.equal(final.reference.deltaX, preview.reference.deltaX);
  assert.notEqual(final.signature, preview.signature);
});
