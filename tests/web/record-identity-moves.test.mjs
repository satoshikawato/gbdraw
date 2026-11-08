// OV-1949 (record-identity-20261027 DESIGN 2-3, A1, A3): a record's identity
// does not depend on how it is drawn, and each Generate moves the rows of every
// record whose request key changes. These are the state-free answers of
// services/record-draw-selection.js.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  drawnRecordIdentities,
  expandedRecordKey,
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
  assert.deepEqual([...legacyCircularRecordKeys(sources)].map(([identity, keys]) => [identity, keys.filter((key) => key.startsWith('record-'))]),
    [['record-1', []], ['record-2:1', ['record-2']], ['record-2:2', ['record-3']]]);
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

// A Circular record of an older Session (50 and older) may have rows under the
// keys of then: `record-<n>` and `circular-<record ID>-#<n>` (`legacyCircularRecordKeys`).
const OLDER = circular({ K: ['A', 'B'] });
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
