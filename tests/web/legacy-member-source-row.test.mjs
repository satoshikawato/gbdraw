// OV-3539: before the first Generate of a Session before 40, the Similarity group members resolve the features Load
// read again from the sources. Load gives each member its canonical identity (the committed request's record key at
// its record index, and its stable ID as its biological feature ID: `migrateLegacyOrthogroupMembers`), while the
// features read again keep the names an older Session used (record index and stable ID) and their own edit keys. The
// member's source row is found by one rule: the same canonical names, read for the feature from the committed
// request, so the drawer and every member reader read the source strand and location, not the member's saved copy,
// which a Session 39 saved in the drawn direction (OV-3442).
import assert from 'node:assert/strict';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };

const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
const { materializeLegacySimilarityAlignment, migrateLegacyOrthogroupMembers } = await import(
  '../../gbdraw/web/js/services/legacy-similarity-alignment.js'
);
const { featureOverrideKey } = await import('../../gbdraw/web/js/services/feature-override-identity.js');
const { ref, reactive } = Vue;

const records = [{ recordKey: 'linear-seq-a' }, { recordKey: 'linear-seq-b' }];
// The features Load read again: each names its drawn record (record index 1 is the second card) and its source
// hash; the second record is drawn reverse complemented, so its features carry another rendered ID.
const sourceRow = (stable, recordIndex, strand, start, end, extra = {}) => ({
  id: `file${recordIndex}_f0`, svg_id: stable, stable_svg_id: stable, stable_feature_id: stable,
  fileIdx: recordIndex, record_idx: recordIndex, feature_index: 12, record_id: `REC${recordIndex}`,
  type: 'CDS', start, end, strand, gene: 'racM', ...extra
});
const livA = sourceRow('stable-liva', 0, '-', 100, 400);
const racM = sourceRow('stable-racm', 1, '-', 3520, 5023, { rendered_svg_id: 'drawn-racm_record_2', renderedSvgId: 'drawn-racm_record_2' });
// The members as a Session 39 saved them: record index, source feature index, stable ID, and the strand in the
// drawn direction.
const savedMember = (stable, recordIndex, strand, start, end) => ({
  recordIndex, featureIndex: 52, featureSvgId: stable, stableFeatureSvgId: stable, start, end, strand, gene: 'racM'
});
const [group] = migrateLegacyOrthogroupMembers([{ id: 'og_18', members: [
  savedMember('stable-liva', 0, '-', 100, 400),
  savedMember('stable-racm', 1, '+', 3520, 5023)
] }], records);

const editorFor = (biologicalFeatures, request = { records }) => {
  const state = {
    mode: 'linear', orthogroups: ref([group]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('og_18'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    // Load installs the features read again as both lists; a drawn one names its drawn element in `svg_id`.
    svgContainer: ref(null), linearSeqs: reactive([]),
    extractedFeatures: ref(biologicalFeatures.map((row) => ({ ...row, svg_id: row.rendered_svg_id ?? row.svg_id }))),
    biologicalFeatures: ref(biologicalFeatures)
  };
  const getCommittedRequest = typeof request === 'function' ? request : () => request;
  return createOrthogroupEditor({ state: withDrawings(state), getCommittedRequest, memberAction: () => null });
};

test('Load gives the saved members their canonical identity from the committed request', () => {
  assert.deepEqual(group.members.map((member) => [member.recordKey, member.biologicalFeatureId]),
    [['linear-seq-a', 'stable-liva'], ['linear-seq-b', 'stable-racm']]);
});

test('a member of a Session before 40 resolves the feature read again, by its canonical identity', () => {
  const members = editorFor([livA, racM]).getEnrichedOrthogroupMembers(group);
  assert.deepEqual(members.map((member) => member.sequenceFeature ?? null), [livA, racM]);
  // The source strand, not the saved drawn one.
  assert.equal(members[1].sequenceFeature?.strand, '-');
  // The features keep their own names: their edit keys stay those the Load migrated the edits to.
  assert.equal(Object.hasOwn(racM, 'recordKey'), false);
  assert.equal(featureOverrideKey(racM), 'file1_f0');
});

test('a feature whose record index names no record of the committed request resolves no member', () => {
  const members = editorFor([livA, racM], { records: records.slice(0, 1) }).getEnrichedOrthogroupMembers(group);
  assert.deepEqual(members.map((member) => member.sequenceFeature ?? null), [livA, null]);
});

// The saved alignment of such a Session names its members by the same rule: a member that saved its stable ID under
// any of its keys gets the names Load gives it, and its source feature index must still be saved.
test('a saved alignment reads each member by the rule Load gives the members their names', () => {
  for (const stableKey of ['stableFeatureSvgId', 'featureSvgId', 'stableFeatureId', 'stable_feature_id']) {
    const members = [
      { recordIndex: 0, featureIndex: 3, [stableKey]: 'stable-liva', representative: true },
      { recordIndex: 1, featureIndex: 7, [stableKey]: 'stable-racm' }
    ];
    const plan = materializeLegacySimilarityAlignment({
      target: 'og_18', records, featureCatalog: null, legacyOrthogroupState: { groups: [{ id: 'og_18', members }] }
    });
    const [migrated] = migrateLegacyOrthogroupMembers([{ id: 'og_18', members }], records);
    assert.deepEqual(plan.records.map(({ anchor }) => [anchor.recordKey, anchor.biologicalFeatureId]),
      migrated.members.map((member) => [member.recordKey, member.biologicalFeatureId]), stableKey);
    assert.deepEqual(plan.records.map(({ status }) => status), ['reference', 'aligned'], stableKey);
  }
  const withoutIndex = [{ recordIndex: 0, stableFeatureId: 'stable-liva', representative: true }];
  assert.throws(() => materializeLegacySimilarityAlignment({
    target: 'og_18', records, featureCatalog: null, legacyOrthogroupState: { groups: [{ id: 'og_18', members: withoutIndex }] }
  }), /lacks a source feature index/);
  const outOfRange = [{ recordIndex: 2, featureIndex: 3, stableFeatureId: 'stable-liva', representative: true }];
  assert.throws(() => materializeLegacySimilarityAlignment({
    target: 'og_18', records, featureCatalog: null, legacyOrthogroupState: { groups: [{ id: 'og_18', members: outOfRange }] }
  }), /no stable record mapping/);
});

// The lookup reads the committed request's record keys only: a committed change that keeps them (a placement, a
// setting) does not read every row again; one that changes them does.
test('the row lookup is rebuilt when the record keys change, not on every committed change', () => {
  let reads = 0;
  const counted = { ...racM };
  delete counted.stable_feature_id;
  Object.defineProperty(counted, 'stable_feature_id', {
    enumerable: true, get: () => { reads += 1; return 'stable-racm'; }
  });
  const committed = ref({ records, layout: { recordGapPx: 30 } });
  const editor = editorFor([livA, counted], () => committed.value);
  const strands = () => editor.getEnrichedOrthogroupMembers(group).map((member) => member.sequenceFeature?.strand);
  assert.deepEqual(strands(), ['-', '-']);
  const built = reads;
  committed.value = { records: records.map((record) => ({ ...record })), layout: { recordGapPx: 40 } };
  assert.deepEqual(strands(), ['-', '-']);
  assert.equal(reads, built, 'same record keys: the rows are not read again');
  committed.value = { records: records.slice(0, 1) };
  assert.deepEqual(strands(), ['-', undefined]);
  assert.ok(reads > built, 'other record keys: the rows are read again');
});
