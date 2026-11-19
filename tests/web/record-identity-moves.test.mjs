// OV-1949 (record-identity-20261027 DESIGN 2-3, A1, A3): a record's identity
// does not depend on how it is drawn, and each Generate moves the rows of every
// record whose request key changes. These are the state-free answers of
// services/record-draw-selection.js.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  cardsTakeLegacyRecordPositions,
  drawnRecordIdentities,
  expandedRecordKey,
  heldCardPrefix,
  legacyCircularRecordKeys,
  recordIdentity,
  recordKeyExpansion,
  recordKeyMoves
} from '../../gbdraw/web/js/services/record-draw-selection.js';
import { buildDisambiguatedRecordEntries } from '../../gbdraw/web/js/services/record-options.js';
import { moveInstanceRuleRecords, namedFeatureInstance as named } from '../../gbdraw/web/js/services/feature-selector.js';
import { movedIdentityValues, movedRecordEdits } from '../../gbdraw/web/js/services/feature-visibility.js';
import { requestFeatureOverrides, requestFeaturePlacements } from '../../gbdraw/web/js/services/feature-placement.js';
import { annotationOptionsPayload } from '../../gbdraw/web/js/services/annotation-state.js';

const entries = (...recordIds) => buildDisambiguatedRecordEntries(recordIds.map((recordId) => ({ recordId })));
const cardsOf = (cards) => new Map(Object.entries(cards).map(([uid, ids]) => [uid, entries(...ids)]));
const circular = (cards) => ({ mode: 'circular', cards: cardsOf(cards) });
const linear = (cards) => ({ mode: 'linear', cards: cardsOf(cards) });
const moves = (previous, next) => Object.fromEntries(recordKeyMoves(previous, next));
const pairs = (drawn) => drawn.map(({ requestKey, identity }) => [requestKey, identity]);

test('a record identity names the record by its input and its position in the file (D-02, D-04)', () => {
  // One rule for both modes: a Circular source card is a card like a Linear one.
  assert.equal(recordIdentity({ uid: 'C', sourceIndex: 1, recordCount: 3 }), 'C:2');
  assert.equal(recordIdentity({ uid: 'U', sourceIndex: 1, recordCount: 2 }), 'U:2');
  // A file of one record keeps the card uid, the key 0.14.2 wrote.
  assert.equal(recordIdentity({ uid: 'U', sourceIndex: 0, recordCount: 1 }), 'U');
  assert.equal(expandedRecordKey('record-1', 2), 'record-1:2');
  assert.equal(expandedRecordKey('record-1', null), 'record-1');
  assert.equal(recordKeyExpansion('U:12', 'U'), 12);
  assert.equal(recordKeyExpansion('U:0', 'U'), null);
  assert.equal(recordKeyExpansion('U:x', 'U'), null);
  assert.equal(recordKeyExpansion('UU:1', 'U'), null);
});

test('Circular: the canvas, one record drawn alone, and a region name each record by its identity', () => {
  const sources = circular({ C: ['A', 'B', 'C'] });
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'C', cardinality: 'all', windows: [] }], sources)),
    [['C:1', 'C:1'], ['C:2', 'C:2'], ['C:3', 'C:3']]);
  const single = [{ recordKey: 'C:2', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'B' }, windows: [] }];
  assert.deepEqual(pairs(drawnRecordIdentities(single, sources)), [['C:2', 'C:2']]);
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'C', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 2 },
    windows: [{ start: 1, end: 10, reverseComplement: false, column: null }] }], sources)),
    [['C', 'C:3']]);
  // A record ID that two records share names none; an unknown card names none.
  assert.deepEqual(drawnRecordIdentities([{ recordKey: 'C', selector: { kind: 'recordId', value: 'A' }, windows: [] }], circular({ C: ['A', 'A'] })), []);
  assert.deepEqual(drawnRecordIdentities(single, circular({})), []);
});

test('Circular canvas <-> one record moves no record row; the keys of an older Session never take a live key', () => {
  const sources = circular({ 'record-1': ['A'], 'record-2': ['B', 'C'] });
  const canvas = drawnRecordIdentities([{ recordKey: 'record-1', cardinality: 'all', windows: [] }, { recordKey: 'record-2', cardinality: 'all', windows: [] }], sources);
  const single = drawnRecordIdentities([{ recordKey: 'record-2:2', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 }, windows: [] }], sources);
  const own = ['record-1', 'record-2:1', 'record-2:2'];
  for (const [previous, next] of [[canvas, single], [single, canvas], [canvas, canvas]]) {
    assert.deepEqual([...recordKeyMoves(previous, next)].filter(([from, to]) => from !== to && own.includes(from)), []);
  }
  // Session 50 numbered the records over the drawing (`record-<n>`): `record-1` is a live identity here, so only
  // `record-2` and `record-3` are aliases, of the records at those places.
  const loaded = { ...sources, legacyOrder: ['record-1', 'record-2'] };
  assert.deepEqual([...legacyCircularRecordKeys(loaded)].map(([identity, keys]) => [identity, keys.filter((key) => key.startsWith('record-'))]),
    [['record-1', []], ['record-2:1', ['record-2']], ['record-2:2', ['record-3']]]);
  // Without a Load-time card order (a Session 51, CLI per-file keys) there is no alias (OV-2867).
  assert.deepEqual([...legacyCircularRecordKeys(sources)], []);
});

test('Linear ALL <-> one record <-> region, a file of one record, and a record 2 -> 3 change', () => {
  const sources = linear({ U: ['r1', 'r2', 'r3'], S: ['only'] });
  const all = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', selector: null, windows: [] }], sources);
  assert.deepEqual(pairs(all), [['U:1', 'U:1'], ['U:2', 'U:2'], ['U:3', 'U:3']]);
  const one = (index) => drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index } }], sources);
  assert.deepEqual(pairs(one(1)), [['U', 'U:2']]);
  assert.deepEqual(moves(all, one(1)), { 'U:2': 'U' });
  // Record 2 -> record 3: both moves apply together.
  assert.deepEqual(moves(one(1), one(2)), { U: 'U:2', 'U:3': 'U' });
  // Record 2 -> record 1, then back to every record.
  assert.deepEqual(moves(one(1), one(0)), { U: 'U:2', 'U:1': 'U' });
  assert.deepEqual(moves(one(0), all), { U: 'U:1' });
  // A cropped record draws the record its selector names.
  const region = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r2' }, windows: [{ start: 1, end: 10, reverseComplement: false, column: null }] }], sources);
  assert.deepEqual(pairs(region), [['U', 'U:2']]);
  // A file of one record has one key in every view.
  const single = drawnRecordIdentities([{ recordKey: 'S', cardinality: 'all', windows: [] }], sources);
  assert.deepEqual(pairs(single), [['S', 'S']]);
  assert.deepEqual(moves(single, drawnRecordIdentities([{ recordKey: 'S', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 0 } }], sources)), {});
  // A split ALL input keys each record `<uid>:<k>`: the identity.
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'U:2', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } }], sources)), [['U:2', 'U:2']]);
  // A removed card (no list) has no identity: its rows keep their key.
  assert.deepEqual(drawnRecordIdentities([{ recordKey: 'gone', cardinality: 'all', windows: [] }], sources), []);
});

// R2-H1: the cards of one uploaded file (U record 1, X record 2) are two
// cards; one card's rows never move onto the other card, also while the
// first card's record 2 waits under its identity.
test('Linear: the rows of one card of a file never move onto another card of that file', () => {
  // A split origin (removed) named U's record 2 as X's; the sources may carry it.
  const sources = { ...linear({ U: ['r1', 'r2'], X: ['r1', 'r2'] }), linearOrigins: new Map([['X', { uid: 'U', sourceIndex: 1 }]]) };
  const draw = (uRecord) => drawnRecordIdentities([
    { recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordId', value: uRecord } },
    { recordKey: 'X', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r2' } }
  ], sources);
  // U on record 2 (its feature edited), Generate; U back on record 1, Generate; Generate again.
  const onRecord2 = draw('r2');
  const backOnRecord1 = draw('r1');
  assert.deepEqual(moves(onRecord2, backOnRecord1), { U: 'U:2', 'U:1': 'U', X: 'X', 'X:2': 'X' });
  assert.equal(recordKeyMoves(backOnRecord1, backOnRecord1).get('U:2'), undefined);
  [moves(onRecord2, backOnRecord1), moves(backOnRecord1, backOnRecord1)].forEach((step) => {
    assert.equal(Object.entries(step).some(([from, to]) => from.startsWith('U') && to === 'X'), false);
  });
});

test('duplicates stay separate: the same file in two cards, and records with one ID', () => {
  const sources = linear({ U: ['r1', 'r2'], V: ['r1', 'r2'] });
  const previous = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', windows: [] }, { recordKey: 'V', cardinality: 'all', windows: [] }], sources);
  const next = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } },
    { recordKey: 'V', cardinality: 'all', windows: [] }], sources);
  // Only card U's record 2 moves; card V's rows stay under V's keys.
  assert.deepEqual(moves(previous, next), { 'U:2': 'U' });
  // Two records with one ID are told apart by position.
  const twins = linear({ W: ['same', 'same'] });
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'W', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } }], twins)), [['W', 'W:2']]);
});

// #494: a replaced source no longer binds its records; the prune decides
// their rows, so no row leaves or enters a key whose source the Generate replaces.
test('no row moves out of or into a record key whose source the Generate replaces', () => {
  // A card keeps its uid: a row moves neither out of its replaced key nor into it.
  const sources = linear({ U: ['r1', 'r2'] });
  const all = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', windows: [] }], sources);
  const one = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } }], sources);
  assert.deepEqual(moves(all, one), { 'U:2': 'U' });
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all, one, new Set(['U:2']))), {});
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all, one, new Set(['U']))), {});
});

// OV-2882 (D-20): a Replace in the same card keeps its first record. The record
// of a one-record file is `U`, the first record of a file of several `U:1`;
// across the Replace they are one record, so its rows move from `U` to `U:1`
// (and back), the replaced key included: the prune then keeps the edits of the
// features the new file still has.
test('the record of a one-record file and the first record of its replacement are one record', () => {
  const one = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', windows: [] }, { recordKey: 'V', cardinality: 'all', windows: [] }],
    circular({ U: ['r1'], V: ['v'] }));
  const several = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', windows: [] }, { recordKey: 'V', cardinality: 'all', windows: [] }],
    circular({ U: ['r1', 'r2'], V: ['v'] }));
  assert.deepEqual(Object.fromEntries(recordKeyMoves(one, several, new Set(['U']))), { U: 'U:1' });
  assert.deepEqual(Object.fromEntries(recordKeyMoves(several, one, new Set(['U:1', 'U:2']))), { 'U:1': 'U' });
  // Another card's uid is never taken for it.
  const other = drawnRecordIdentities([{ recordKey: 'W', cardinality: 'all', windows: [] }], circular({ W: ['w1', 'w2'] }));
  assert.deepEqual(Object.fromEntries(recordKeyMoves(one, other)), {});
});

// OV-2867: an older Session's flat positions count over its cards in their
// order at Load, never in the order of now: a card moved after Load keeps its
// records' aliases, and the walk stops at the first card gone or whose file
// changed, so a removed card retargets no alias.
test('legacy aliases count over the cards as loaded and stop at the first one gone or replaced', () => {
  const fileA = Object.freeze({ fileId: 'file-a' });
  const fileB = Object.freeze({ fileId: 'file-b' });
  const marks = [{ uid: 'source-1', primaryFile: fileA, pairedFile: null, recordIds: null },
    { uid: 'source-2', primaryFile: fileB, pairedFile: null, recordIds: null }];
  const input = (uid, primaryFile) => ({ uid, primaryFile, pairedFile: null });
  // File 2 moved up: the order of Load stays.
  assert.deepEqual(heldCardPrefix(marks, [input('source-2', fileB), input('source-1', fileA)]), ['source-1', 'source-2']);
  // File 1 removed, or replaced: no card of then is placed past it.
  assert.deepEqual(heldCardPrefix(marks, [input('source-2', fileB)]), []);
  assert.deepEqual(heldCardPrefix(marks, [input('source-1', Object.freeze({ fileId: 'file-c' })), input('source-2', fileB)]), []);
  assert.deepEqual(heldCardPrefix(marks, [input('source-1', fileA)]), ['source-1']);
  const moved = { ...circular({ 'source-2': ['b'], 'source-1': ['a'] }), legacyOrder: ['source-1', 'source-2'] };
  assert.deepEqual(Object.fromEntries([...legacyCircularRecordKeys(moved)].map(([identity, keys]) => [identity, keys[0]])),
    { 'source-1': 'record-1', 'source-2': 'record-2' });
  // A card not read yet ends the walk.
  assert.deepEqual([...legacyCircularRecordKeys({ ...circular({ 'source-1': [], 'source-2': ['b'] }), legacyOrder: ['source-1', 'source-2'] })], []);
  // Which migrated cards take flat positions (D-12): one card, or several not keyed by the CLI's per-file keys.
  assert.equal(cardsTakeLegacyRecordPositions([{ uid: 'record-1' }]), true);
  assert.equal(cardsTakeLegacyRecordPositions([{ uid: 'source-1' }, { uid: 'source-2' }]), true);
  assert.equal(cardsTakeLegacyRecordPositions([{ uid: 'record-1' }, { uid: 'record-2' }]), false);
  assert.equal(cardsTakeLegacyRecordPositions([]), false);
});

// A Circular record of an older Session (50 and older) may have rows under the
// keys of then: `record-<n>` and `circular-<record ID>-#<n>` (`legacyCircularRecordKeys`).
const OLDER = { ...circular({ K: ['A', 'B'] }), legacyOrder: ['K'] };
const recordB = drawnRecordIdentities([{ recordKey: 'K:2', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'B' } }], OLDER);

test('aliases bring the rows an older Session left under its keys, after the record\'s own rows', () => {
  assert.deepEqual([...recordKeyMoves(recordB, recordB)], [['K:2', 'K:2'], ['record-2', 'K:2'], ['circular-B-#2', 'K:2']]);
});

// OV-1949: a Generate that moves a record's rows renames the record key of each
// `hash` rule naming one of its feature instances, keeping how the rule was
// written; Python compares an instance without letter case.
test('moveInstanceRuleRecords renames the instance rules of moved record keys only', () => {
  const rules = [
    { qual: 'hash', val: 'f1@record-2', color: '#f00' },
    { qual: 'hash', val: '^f2@record-2$', color: '#0f0' },
    { qual: 'hash', val: 'f3@RECORD-2', color: '#00f' },
    { qual: 'hash', val: 'f4@record-1', color: '#000' },
    { qual: 'hash', val: '^f5@record-.*$', color: '#111' },
    { qual: 'product', val: 'f6@record-2', color: '#222' }
  ];
  const moved = moveInstanceRuleRecords(rules, new Map([['record-2', 'circular-B-#2']]));
  assert.deepEqual(moved.map((rule) => rule.val), [
    'f1@circular-B-#2', '^f2@circular-B-#2$', 'f3@circular-B-#2', 'f4@record-1', '^f5@record-.*$', 'f6@record-2'
  ]);
  assert.equal(named(moved[1])?.recordKey, 'circular-B-#2');
  // Nothing moved: the same list.
  assert.equal(moveInstanceRuleRecords(rules, new Map()), rules);
  assert.equal(moveInstanceRuleRecords(rules, new Map([['other', 'x']])), rules);
  // The composer writes an exact pattern on request.
  assert.equal(named(rules[0])?.renamed('record-1:2', 'exact'), '^f1@record-1:2$');
});

test('two rules that name one feature instance after a move keep the rule of the record\'s own key', () => {
  const rules = [
    { qual: 'hash', val: '^f@record-2$', color: '#f00' },
    { qual: 'hash', val: '^f@K:2$', color: '#00f' },
    { qual: 'hash', val: '^g@record-2$', color: '#0f0' },
    { qual: 'hash', val: '^h@K:1$', color: '#111' },
    { qual: 'hash', val: '^h@K:1$', color: '#222' }
  ];
  // The record's own key wins: blue stays, red goes; two rules that already named one instance both stay.
  assert.deepEqual(moveInstanceRuleRecords(rules, recordKeyMoves(recordB, recordB)).map((rule) => [rule.val, rule.color]), [
    ['^f@K:2$', '#00f'], ['^g@K:2$', '#0f0'], ['^h@K:1$', '#111'], ['^h@K:1$', '#222']
  ]);
});

test('a row under a record\'s own key wins over an alias row field by field', () => {
  const key = (recordKey) => JSON.stringify([recordKey, 'f']);
  const row = (recordKey, fields) => ({ recordKey, biologicalFeatureId: 'f', featureVisibility: null,
    labelVisibility: null, labelText: null, labelSourceText: null, ...fields });
  const edits = {
    featureOverrides: {
      [key('K:2')]: row('K:2', { labelText: 'own' }),
      [key('record-2')]: row('record-2', { labelText: 'alias', featureVisibility: 'off' })
    },
    featurePlacementOverrides: {}, featureStrokeOverrides: {}, featureColorOverrides: {}, annotationSets: []
  };
  assert.deepEqual(movedRecordEdits(edits, recordKeyMoves(recordB, recordB)).featureOverrides,
    { [key('K:2')]: row('K:2', { labelText: 'own', featureVisibility: 'off' }) });
});

test('the request of a view-switch Generate reads the moved rows; the drawing keeps its own until commit', () => {
  const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
  const drafts = {
    featureOverrides: { [key('record-2', 'f1')]: { recordKey: 'record-2', biologicalFeatureId: 'f1', featureVisibility: 'off',
      labelVisibility: null, labelText: 'edited', labelSourceText: null } },
    featurePlacementOverrides: { [key('record-2', 'f1')]: { recordKey: 'record-2', biologicalFeatureId: 'f1', placement: { kind: 'main' } } },
    featureStrokeOverrides: { 'record-2\0f1': { width: 2 } },
    featureColorOverrides: { 'record-2\0f1': { color: '#ff0000', caption: '' } },
    annotationSets: [{ id: 'set', annotations: [{ target: { kind: 'featureIdentity', recordKey: 'record-2', biologicalFeatureId: 'f1',
      envelope: 'outer_bounds', circularPath: 'shortest' } }] }]
  };
  const before = structuredClone(drafts);
  const moves = new Map([['record-2', 'circular-B-#2']]);
  const moved = movedRecordEdits(drafts, moves);
  const records = [{ recordKey: 'circular-B-#2', cardinality: 'exactly_one' }];
  assert.deepEqual(requestFeatureOverrides(moved.featureOverrides, records).map((row) => [row.recordKey, row.labelText]),
    [['circular-B-#2', 'edited']]);
  assert.deepEqual(requestFeaturePlacements(moved.featurePlacementOverrides, 'circular', records).map((row) => row.recordKey), ['circular-B-#2']);
  assert.deepEqual(Object.keys(moved.featureStrokeOverrides), ['circular-B-#2\0f1']);
  assert.deepEqual(Object.keys(moved.featureColorOverrides), ['circular-B-#2\0f1']);
  assert.deepEqual(annotationOptionsPayload(moved.annotationSets, records).sets[0].annotations.map((item) => item.target.recordKey),
    ['circular-B-#2']);
  // The bulk label projection's keys move the same way.
  assert.deepEqual(movedIdentityValues({ [key('record-2', 'f1')]: 'bulk' }, moves), { [key('circular-B-#2', 'f1')]: 'bulk' });
  // The drawing's own rows are untouched.
  assert.deepEqual(drafts, before);
  // Without the moves the request of the new view has none of the record's rows (the 0.14.2 behavior).
  assert.deepEqual(requestFeatureOverrides(drafts.featureOverrides, records), []);
});

// OV-3846 (D-05, refining D-20; D-23): a Replace in the same card continues
// each record by its accession (its record ID without the version) or its
// record ID when every ID of both files is unique, else by place; a record
// the new file has no such record for leaves. Every row a record holds moves
// by this one rule: the feature edits here (`recordKeyMoves`), the Depth rows
// in `resliceCardRows`.
test('a reordered Replace moves each record\'s rows to the record with its ID; duplicated IDs follow places', async () => {
  const { continuedRecordPlaces } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  assert.deepEqual(continuedRecordPlaces(['A', 'B'], ['B', 'A']), [1, 0]);
  assert.deepEqual(continuedRecordPlaces(['A', 'B'], ['B']), [-1, 0], 'A leaves; B keeps its rows');
  assert.deepEqual(continuedRecordPlaces(['NC_1.1'], ['NC_1.2']), [0], 'a new version of the record keeps it (D-32)');
  assert.deepEqual(continuedRecordPlaces(['A.1', 'B.1'], ['B.2', 'A.1']), [1, 0]);
  assert.deepEqual(continuedRecordPlaces(['A'], ['B']), [-1], 'a record with another ID does not continue A (D-23)');
  assert.deepEqual(continuedRecordPlaces(['A', 'B'], ['A', 'C']), [0, -1]);
  // Two versions of one accession in a file: the exact ID decides.
  assert.deepEqual(continuedRecordPlaces(['A.1'], ['A.2', 'A.1']), [1]);
  assert.deepEqual(continuedRecordPlaces(['A.1', 'A.2'], ['A.3']), [-1, -1]);
  assert.deepEqual(continuedRecordPlaces(['A', 'B'], ['B', 'B']), [0, 1], 'a duplicated ID: by place');
  assert.deepEqual(continuedRecordPlaces(['A', 'A', 'B'], ['B']), [0, -1, -1]);
  const all = (cards) => drawnRecordIdentities(Object.keys(cards).map((recordKey) => ({ recordKey, cardinality: 'all', windows: [] })), circular(cards));
  const replaced = new Set(['U:1', 'U:2']);
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all({ U: ['A', 'B'], V: ['v'] }), all({ U: ['B', 'A'], V: ['v'] }), replaced)),
    { 'U:1': 'U:2', 'U:2': 'U:1' });
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all({ U: ['A', 'B'] }), all({ U: ['B'] }), replaced)), { 'U:2': 'U' });
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all({ U: ['A', 'A'] }), all({ U: ['A', 'A'] }), replaced)), {});
  assert.deepEqual(Object.fromEntries(recordKeyMoves(all({ U: ['A', 'B'] }), all({ U: ['B', 'B'] }), replaced)), {});
  // The edits move together and the prune keeps each one whose feature its record still has.
  const { pruneUnmatchedFeatureOverrides } = await import('../../gbdraw/web/js/services/feature-visibility.js');
  const { FEATURE_EDITS_SLICE } = await import('../../gbdraw/web/js/services/drawing/feature-edits.js');
  const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
  const row = (recordKey, biologicalFeatureId, labelText) => ({
    recordKey, biologicalFeatureId, featureVisibility: null, labelVisibility: null, labelText, labelSourceText: null
  });
  const edits = FEATURE_EDITS_SLICE.create({ featureOverrides: {
    [key('U:1', 'a1')]: row('U:1', 'a1', 'on A'), [key('U:2', 'b1')]: row('U:2', 'b1', 'on B')
  }, featurePlacementOverrides: {} });
  pruneUnmatchedFeatureOverrides({
    featureEdits: edits.handle,
    replacedRecordKeys: [...replaced],
    previousRecords: [{ recordKey: 'U', cardinality: 'all' }],
    currentRecords: [{ recordKey: 'U', cardinality: 'all' }],
    biologicalFeatures: [{ record_key: 'U:1', biological_feature_id: 'b1' }, { record_key: 'U:2', biological_feature_id: 'a1' }],
    moves: recordKeyMoves(all({ U: ['A', 'B'] }), all({ U: ['B', 'A'] }), replaced)
  });
  assert.deepEqual(edits.view.featureOverrides, {
    [key('U:2', 'a1')]: row('U:2', 'a1', 'on A'), [key('U:1', 'b1')]: row('U:1', 'b1', 'on B')
  });
});

// OV-3846 review (MEDIUM 1): several Replaces of a card before one Generate. The Depth rows move at the read after
// each Replace and the feature edits at the Generate; both continue from the records the committed request drew
// (`continuedSinceCommit` for the rows, `recordKeyMoves` for the edits), so the rows of a record never part from its
// edits. Where IDs name records the one-step rule composes; a file between them whose IDs repeat goes by place, and
// there the step-by-step reading would put A's Depth on B while A's edit goes to A.
test('Replaces before one Generate move a record\'s Depth rows where that Generate moves its feature edits', async () => {
  const { continuedRecordPlaces, continuedSinceCommit, resliceCardRows, cardRowsLeaving } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const all = (cards) => drawnRecordIdentities(Object.keys(cards).map((recordKey) => ({ recordKey, cardinality: 'all', windows: [] })), circular(cards));
  const committed = new Map([['U', ['A', 'B']]]);
  const card = (...recordIds) => [{ uid: 'U', recordIds }];
  // (A, B) committed, then (A, A), then (B, A): the rows of (A, A) stood for A and B by place.
  assert.deepEqual(continuedRecordPlaces(['A', 'A'], ['B', 'A']), [0, 1]);
  assert.deepEqual(continuedSinceCommit(['A', 'B'], ['A', 'A'], ['B', 'A']), [1, 0]);
  assert.deepEqual(resliceCardRows([['a.tsv'], ['b.tsv']], card('A', 'A'), card('B', 'A'), committed, () => []), [['b.tsv'], ['a.tsv']]);
  assert.deepEqual(moves(all({ U: ['A', 'B'] }), all({ U: ['B', 'A'] })), { 'U:1': 'U:2', 'U:2': 'U:1' }, 'A\'s edits go to U:2, where A\'s Depth went');
  assert.deepEqual(cardRowsLeaving([['a.tsv'], ['b.tsv']], card('A', 'A'), card('B', 'A'), committed), []);
  // (A.1, B) -> (B, A.2) -> (A.3, B): A's rows go to its newest version, as A's edits do.
  assert.deepEqual(continuedSinceCommit(['A.1', 'B'], ['B', 'A.2'], ['A.3', 'B']), [1, 0]);
  assert.deepEqual(continuedRecordPlaces(['A.1', 'B'], ['A.3', 'B']), [0, 1]);
  // (A, B) -> (B) -> (A, B): B's row stays with B (A's row left at the first Replace, which asked).
  assert.deepEqual(continuedSinceCommit(['A', 'B'], ['B'], ['A', 'B']), [1]);
  // A record the committed request did not draw, or a card it did not draw, continues step by step.
  assert.deepEqual(continuedSinceCommit(['A', 'B'], ['C'], ['B', 'C']), [1], 'C, new since the Generate, keeps its row');
  assert.deepEqual(continuedSinceCommit(null, ['C', 'D'], ['B', 'A']), continuedRecordPlaces(['C', 'D'], ['B', 'A']));
  assert.deepEqual(resliceCardRows([['c.tsv']], [{ uid: 'V', recordIds: ['C.1'] }], [{ uid: 'V', recordIds: ['C.2'] }], committed, () => []), [['c.tsv']]);
  // One record of the new file continues one record: a claim by its exact ID before one by place.
  assert.deepEqual(continuedSinceCommit(['A', 'B'], ['A', 'A'], ['A', 'C']), [0, -1]);
  // With no Replace since the Generate it is the one-step rule.
  for (const [then, next] of [[['A', 'B'], ['B', 'A']], [['A', 'B'], ['B']], [['A.1'], ['A.2']], [['A'], ['B']], [['A', 'B'], ['A', 'C']], [['A', 'B'], ['B', 'B']], [['A', 'A', 'B'], ['B']]]) {
    assert.deepEqual(continuedSinceCommit(then, then, next), continuedRecordPlaces(then, next), `${then} -> ${next}`);
  }
});

// OV-3846 review of 91340c36c (HIGH 1, D-23): two Replaces of one card before a Generate. Each record's Depth row
// (moved at the read after each Replace, after the question) and its feature edits (moved at the Generate) end on the
// same record, and never on a record with another ID: a record continues as the record of the new file with its ID,
// or with its accession without the version; a record the new file has none for leaves, and the question names it.
// The one parting left is the user's answer: a Depth row removed at the question while the record comes back with its
// edits.
test('two Replaces before a Generate never draw a record\'s Depth or edits on a record with another ID', async () => {
  const { cardRowsLeaving, resliceCardRows } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const all = (ids) => drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', windows: [] }], circular({ U: ids }));
  const card = (recordIds) => [{ uid: 'U', recordIds }];
  /** Depth rows named by their record, the records each question names, and the edit moves of the Generate. */
  const replaceTwice = (committedIds, ...files) => {
    const committed = new Map([['U', committedIds]]);
    let rows = committedIds.map((id) => [`${id}.tsv`]);
    let marked = committedIds;
    const asked = [];
    for (const next of files) {
      // The question names the rows that hold a file (`circularDepthFilesLeaving`).
      asked.push(cardRowsLeaving(rows, card(marked), card(next), committed).flatMap(({ row }) => row));
      rows = resliceCardRows(rows, card(marked), card(next), committed, () => []);
      marked = next;
    }
    return { depth: rows.map((row) => row[0] ?? null), asked, edits: moves(all(committedIds), all(marked)) };
  };
  // (a) (A, B) -> (A) -> (B): B's Depth leaves at the first question, A's at the second; B keeps its edits (U:2 -> U).
  assert.deepEqual(replaceTwice(['A', 'B'], ['A'], ['B']), { depth: [null], asked: [['B.tsv'], ['A.tsv']], edits: { 'U:2': 'U' } });
  // (b) (A) -> (B) -> (C, B): A leaves at the first question; nothing of A reaches B or C.
  assert.deepEqual(replaceTwice(['A'], ['B'], ['C', 'B']), { depth: [null, null], asked: [['A.tsv'], []], edits: {} });
  // (c) (A) -> (A, B) -> (B): A's Depth leaves at the second question, and A's edits do not move onto B.
  assert.deepEqual(replaceTwice(['A'], ['A', 'B'], ['B']), { depth: [null], asked: [[], ['A.tsv']], edits: {} });
  // (d) (A.1) -> (C, D) -> (D, A.2): A.1 leaves at the first question; its edits go to its new version A.2.
  assert.deepEqual(replaceTwice(['A.1'], ['C', 'D'], ['D', 'A.2']), { depth: [null, null], asked: [['A.1.tsv'], []], edits: { U: 'U:2' } });
  // A single record replaced by a record with another ID leaves (asked); a new version of it continues.
  assert.deepEqual(replaceTwice(['A'], ['B']), { depth: [null], asked: [['A.tsv']], edits: {} });
  assert.deepEqual(replaceTwice(['NC_1.1'], ['NC_1.2']), { depth: ['NC_1.1.tsv'], asked: [[]], edits: {} });
  // Over every two-Replace sequence of files of one or two records with IDs A to D: no Depth row or edit lands on a
  // record with another ID, and a Depth row and its record's edits (`continuedRecordPlaces`, which `recordKeyMoves`
  // applies) part only where a question removed the row.
  const { continuedRecordPlaces } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const files = ['A', 'B', 'C', 'D'].flatMap((x) => [[x], ...['A', 'B', 'C', 'D'].filter((y) => y !== x).map((y) => [x, y])]);
  let checked = 0;
  for (const a of files) for (const b of files) for (const c of files) {
    const { depth, asked } = replaceTwice(a, b, c);
    const direct = continuedRecordPlaces(a, c);
    a.forEach((id, place) => {
      const label = `${[a, b, c].map((f) => f.join('')).join(' -> ')}: ${id}`;
      const depthPlace = depth.indexOf(`${id}.tsv`);
      if (direct[place] >= 0) assert.equal(c[direct[place]], id, `${label}: edits on another record`);
      if (depthPlace >= 0) assert.deepEqual([c[depthPlace], direct[place]], [id, depthPlace], `${label}: Depth apart from its edits`);
      else assert.ok(direct[place] < 0 || asked.flat().includes(`${id}.tsv`), `${label}: Depth left without a question`);
      checked += 1;
    });
  }
  assert.equal(checked, 7168);
  // Where a file between them or the new file repeats an ID (by place), a Depth row still goes only where its record's
  // edits go, or leaves through the question.
  const repeated = [...files, ...['A', 'B'].map((x) => [x, x])];
  for (const a of files) for (const b of repeated) for (const c of repeated) {
    const { depth, asked } = replaceTwice(a, b, c);
    const direct = continuedRecordPlaces(a, c);
    a.forEach((id, place) => {
      const label = `${[a, b, c].map((f) => f.join('')).join(' -> ')}: ${id}`;
      const depthPlace = depth.indexOf(`${id}.tsv`);
      if (depthPlace >= 0) assert.equal(direct[place], depthPlace, `${label}: Depth apart from its edits`);
      else assert.ok(direct[place] < 0 || asked.flat().includes(`${id}.tsv`), `${label}: Depth left without a question`);
    });
  }
});

// A Linear card's request record names its record by ID (`region_record_id`) under the card's uid: a reordered
// Replace keeps drawing that record under the same request key, so its rows stay where they are.
test('Linear: a reordered Replace keeps the rows of the record a card draws by its ID', () => {
  const one = (id) => [{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordId', value: id }, windows: [] }];
  const then = drawnRecordIdentities(one('A'), linear({ U: ['A', 'B'] }));
  const now = drawnRecordIdentities(one('A'), linear({ U: ['B', 'A'] }));
  assert.deepEqual(pairs(then), [['U', 'U:1']]);
  assert.deepEqual(pairs(now), [['U', 'U:2']]);
  // The rows under U stay (U still draws A); B, not drawn, takes the rows under its identity along.
  assert.deepEqual(moves(then, now), { 'U:2': 'U:1' });
});
