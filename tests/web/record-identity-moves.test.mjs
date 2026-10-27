// OV-1949 (record-identity-20261027 DESIGN 2-3, A1, A3): a record's identity
// does not depend on how it is drawn, and each Generate moves the rows of every
// record whose request key changes. These are the state-free answers of
// services/record-draw-selection.js.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  circularSingleRecordKey,
  drawnRecordIdentities,
  expandedRecordKey,
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
const circular = (...recordIds) => ({ mode: 'circular', circular: entries(...recordIds) });
const linear = (cards) => ({ mode: 'linear', linear: new Map(Object.entries(cards).map(([uid, ids]) => [uid, entries(...ids)])) });
const moves = (previous, next) => Object.fromEntries(recordKeyMoves(previous, next));
const pairs = (drawn) => drawn.map(({ requestKey, identity }) => [requestKey, identity]);

test('a record identity names the record by its input and its position in the file (D-02, D-04)', () => {
  assert.equal(recordIdentity({ mode: 'circular', uid: '', sourceIndex: 1, recordCount: 3 }), 'record-2');
  assert.equal(recordIdentity({ mode: 'linear', uid: 'U', sourceIndex: 1, recordCount: 2 }), 'U:2');
  // A file of one record keeps the card uid, the key 0.14.2 wrote.
  assert.equal(recordIdentity({ mode: 'linear', uid: 'U', sourceIndex: 0, recordCount: 1 }), 'U');
  assert.equal(expandedRecordKey('record-1', 2), 'record-1:2');
  assert.equal(expandedRecordKey('record-1', null), 'record-1');
  assert.equal(recordKeyExpansion('U:12', 'U'), 12);
  assert.equal(recordKeyExpansion('U:0', 'U'), null);
  assert.equal(recordKeyExpansion('U:x', 'U'), null);
  assert.equal(recordKeyExpansion('UU:1', 'U'), null);
});

test('a Circular record drawn alone does not keep a preserved key that names another record (A3)', () => {
  // A CLI Session keyed its only input `record-1`; that is record 1's identity, not record 3's.
  assert.equal(circularSingleRecordKey({ recordKey: 'record-1', recordId: 'C', selector: '#3', sourceIndex: 2 }), 'circular-C-#3');
  assert.equal(circularSingleRecordKey({ recordKey: 'record-3', recordId: 'C', selector: '#3', sourceIndex: 2 }), 'record-3');
  assert.equal(circularSingleRecordKey({ recordKey: 'python-key', recordId: 'C', selector: '#3', sourceIndex: 2 }), 'python-key');
});

test('Circular: the canvas, one record drawn alone, and a loaded canvas numbered in request order', () => {
  const sources = circular('A', 'B', 'C');
  const canvas = [{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 0 } },
    { recordKey: 'record-2', selector: { kind: 'recordIndex', index: 1 } }];
  assert.deepEqual(pairs(drawnRecordIdentities(canvas, sources)), [['record-1', 'record-1'], ['record-2', 'record-2']]);
  const single = [{ recordKey: 'circular-B-#2', selector: { kind: 'recordId', value: 'B' } }];
  const singleDrawn = drawnRecordIdentities(single, sources);
  assert.deepEqual(pairs(singleDrawn), [['circular-B-#2', 'record-2']]);
  // Its aliases are the two keys 0.14.2 gave it.
  assert.deepEqual(singleDrawn[0].aliases, ['record-2', circularSingleRecordKey(sources.circular[1])]);
  // A region names its record through the region's selector.
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'k', region: { selector: { kind: 'recordIndex', index: 2 } } }], sources)),
    [['k', 'record-3']]);
  // A loaded canvas numbers its keys in request order; the selector names the record.
  const loaded = [{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 2 } },
    { recordKey: 'record-2', selector: { kind: 'recordIndex', index: 0 } }];
  assert.deepEqual(pairs(drawnRecordIdentities(loaded, sources)), [['record-1', 'record-3'], ['record-2', 'record-1']]);
  // A record ID that two records share names none; an unknown list names none.
  assert.deepEqual(drawnRecordIdentities([{ recordKey: 'x', selector: { kind: 'recordId', value: 'A' } }], circular('A', 'A')), []);
  assert.deepEqual(drawnRecordIdentities(single, circular()), []);
});

test('Circular canvas <-> one record <-> batch moves the record rows both ways', () => {
  const sources = circular('A', 'B');
  const canvas = drawnRecordIdentities([{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 0 } },
    { recordKey: 'record-2', selector: { kind: 'recordIndex', index: 1 } }], sources);
  const single = drawnRecordIdentities([{ recordKey: 'circular-B-#2', selector: { kind: 'recordId', value: 'B' } }], sources);
  // Canvas -> record 2 alone: record 2 moves, record 1 already waits under its identity.
  assert.deepEqual(moves(canvas, single), { 'record-2': 'circular-B-#2' });
  // Back to the canvas (a batch has the same keys).
  // Back on the canvas, each record's alias is its key drawn alone, where
  // 0.14.2 may have left rows; a record's rows already under its key come
  // before them (a move to itself).
  assert.deepEqual([...recordKeyMoves(single, canvas)],
    [['circular-B-#2', 'record-2'], ['record-1', 'record-1'], ['circular-A-#1', 'record-1']]);
  assert.deepEqual([...recordKeyMoves(canvas, canvas)], [
    ['record-1', 'record-1'], ['circular-A-#1', 'record-1'], ['record-2', 'record-2'], ['circular-B-#2', 'record-2']
  ]);
});

test('A3: preserved CLI keys of another record are permuted together; aliases never take another record\'s key', () => {
  const sources = circular('A', 'B', 'C');
  // A CLI Session drew R3 as `record-1` and R1 as `record-2`; a Web canvas keys them by file position.
  const loaded = drawnRecordIdentities([{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 2 } },
    { recordKey: 'record-2', selector: { kind: 'recordIndex', index: 0 } }], sources);
  const canvas = drawnRecordIdentities([0, 2].map((index) => ({
    recordKey: `record-${index + 1}`, selector: { kind: 'recordIndex', index }
  })), sources);
  assert.deepEqual(moves(loaded, canvas), {
    'record-1': 'record-3', 'record-2': 'record-1', 'circular-A-#1': 'record-1', 'circular-C-#3': 'record-3'
  });
  // R2 is drawn in neither: its identity `record-2` is R1's previous key, so it does not move there.
  assert.equal(recordKeyMoves(loaded, canvas).get('record-2'), 'record-1');
});

test('Linear ALL <-> one record <-> region, a file of one record, and a record 2 -> 3 change', () => {
  const sources = linear({ U: ['r1', 'r2', 'r3'], S: ['only'] });
  const all = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', selector: null, region: null }], sources);
  assert.deepEqual(pairs(all), [['U:1', 'U:1'], ['U:2', 'U:2'], ['U:3', 'U:3']]);
  const one = (index) => drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index } }], sources);
  assert.deepEqual(pairs(one(1)), [['U', 'U:2']]);
  assert.deepEqual(moves(all, one(1)), { 'U:2': 'U' });
  // Record 2 -> record 3: both moves apply together.
  assert.deepEqual(moves(one(1), one(2)), { U: 'U:2', 'U:3': 'U' });
  // Record 2 -> record 1, then back to every record.
  assert.deepEqual(moves(one(1), one(0)), { U: 'U:2', 'U:1': 'U' });
  assert.deepEqual(moves(one(0), all), { U: 'U:1' });
  // A region draws the record its selector names.
  const region = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', region: { selector: { kind: 'recordId', value: 'r2' } } }], sources);
  assert.deepEqual(pairs(region), [['U', 'U:2']]);
  // A file of one record has one key in every view.
  const single = drawnRecordIdentities([{ recordKey: 'S', cardinality: 'all' }], sources);
  assert.deepEqual(pairs(single), [['S', 'S']]);
  assert.deepEqual(moves(single, drawnRecordIdentities([{ recordKey: 'S', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 0 } }], sources)), {});
  // A split ALL input keys each record `<uid>:<k>`: the identity.
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'U:2', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } }], sources)), [['U:2', 'U:2']]);
  // A removed card (no list) has no identity: its rows keep their key.
  assert.deepEqual(drawnRecordIdentities([{ recordKey: 'gone', cardinality: 'all' }], sources), []);
});

test('Linear: a Generate that splits a card drawing every record moves each record to its own card', () => {
  // A loaded Session drew U's two records as U:1 and U:2; Generate gives
  // record 2 the new card X (`expandDiscoveredLinearRecords`).
  const before = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all', selector: null, region: null }],
    linear({ U: ['r1', 'r2'] }));
  const split = { ...linear({ U: ['r1', 'r2'], X: ['r1', 'r2'] }), linearOrigins: new Map([['X', { uid: 'U', sourceIndex: 1 }]]) };
  const after = drawnRecordIdentities([
    { recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r1' } },
    { recordKey: 'X', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r2' } }
  ], split);
  assert.deepEqual(moves(before, after), { 'U:1': 'U', 'X:2': 'X', 'U:2': 'X' });
  // Without the origin, record 2's rows wait under U:2.
  assert.deepEqual(moves(before, drawnRecordIdentities([
    { recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r1' } },
    { recordKey: 'X', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r2' } }
  ], linear({ U: ['r1', 'r2'], X: ['r1', 'r2'] }))), { 'U:1': 'U', 'X:2': 'X' });
  // The origin names X's record only while X draws that record, and never a
  // key that U draws again.
  const otherRecord = drawnRecordIdentities([{ recordKey: 'X', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r1' } }], split);
  assert.deepEqual(otherRecord[0].aliases, ['X:1']);
  const uAllAgain = drawnRecordIdentities([
    { recordKey: 'U', cardinality: 'all', selector: null, region: null },
    { recordKey: 'X', cardinality: 'exactly_one', selector: { kind: 'recordId', value: 'r2' } }
  ], split);
  assert.equal(recordKeyMoves(after, uAllAgain).get('U:2'), undefined);
});

test('duplicates stay separate: the same file in two cards, and records with one ID', () => {
  const sources = linear({ U: ['r1', 'r2'], V: ['r1', 'r2'] });
  const previous = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'all' }, { recordKey: 'V', cardinality: 'all' }], sources);
  const next = drawnRecordIdentities([{ recordKey: 'U', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } },
    { recordKey: 'V', cardinality: 'all' }], sources);
  // Only card U's record 2 moves; card V's rows stay under V's keys.
  assert.deepEqual(moves(previous, next), { 'U:2': 'U' });
  // Two records with one ID are told apart by position.
  const twins = linear({ W: ['same', 'same'] });
  assert.deepEqual(pairs(drawnRecordIdentities([{ recordKey: 'W', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 1 } }], twins)), [['W', 'W:2']]);
});

test('aliases bring rows a 0.14.2 view switch left dormant, after the record\'s own rows', () => {
  const sources = circular('A', 'B');
  // A 0.14.2 Session saved record 2 drawn alone; its canvas rows still wait under `record-2`.
  const single = drawnRecordIdentities([{ recordKey: 'circular-B-#2', selector: { kind: 'recordId', value: 'B' } }], sources);
  // Generated alone again: the alias `record-2` joins the record's own key,
  // whose rows win field by field.
  assert.deepEqual([...recordKeyMoves(single, single)], [['circular-B-#2', 'circular-B-#2'], ['record-2', 'circular-B-#2']]);
  // Generated on the canvas: the record's own key moves first, the alias is the target itself.
  const canvas = drawnRecordIdentities([{ recordKey: 'record-2', selector: { kind: 'recordIndex', index: 1 } }], sources);
  assert.deepEqual([...recordKeyMoves(single, canvas)], [['circular-B-#2', 'record-2']]);
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
  const sources = circular('A', 'B');
  const single = drawnRecordIdentities([{ recordKey: 'circular-B-#2', selector: { kind: 'recordId', value: 'B' } }], sources);
  const canvas = drawnRecordIdentities([{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 0 } },
    { recordKey: 'record-2', selector: { kind: 'recordIndex', index: 1 } }], sources);
  // A 0.14.2 Session colored f on the canvas (red) and drawn alone (blue).
  const rules = [
    { qual: 'hash', val: '^f@record-2$', color: '#f00' },
    { qual: 'hash', val: '^f@circular-B-#2$', color: '#00f' },
    { qual: 'hash', val: '^g@record-2$', color: '#0f0' },
    { qual: 'hash', val: '^h@record-1$', color: '#111' },
    { qual: 'hash', val: '^h@record-1$', color: '#222' }
  ];
  // On the canvas the canvas key is the record's own: red stays, blue goes;
  // two rules that already named one instance both stay.
  assert.deepEqual(moveInstanceRuleRecords(rules, recordKeyMoves(canvas, canvas)).map((rule) => [rule.val, rule.color]), [
    ['^f@record-2$', '#f00'], ['^g@record-2$', '#0f0'], ['^h@record-1$', '#111'], ['^h@record-1$', '#222']
  ]);
  // Drawn alone again the single key is its own: blue stays, red goes.
  assert.deepEqual(moveInstanceRuleRecords(rules, recordKeyMoves(single, single)).map((rule) => [rule.val, rule.color]), [
    ['^f@circular-B-#2$', '#00f'], ['^g@circular-B-#2$', '#0f0'], ['^h@record-1$', '#111'], ['^h@record-1$', '#222']
  ]);
  // Moved from its key drawn alone to the canvas, blue wins over a row the canvas left.
  assert.deepEqual(moveInstanceRuleRecords(rules, recordKeyMoves(single, canvas)).map((rule) => [rule.val, rule.color]), [
    ['^f@record-2$', '#00f'], ['^g@record-2$', '#0f0'], ['^h@record-1$', '#111'], ['^h@record-1$', '#222']
  ]);
});

test('a row under a record\'s own key wins over an alias row field by field', () => {
  const sources = circular('A', 'B');
  const canvas = drawnRecordIdentities([{ recordKey: 'record-2', selector: { kind: 'recordIndex', index: 1 } }], sources);
  const key = (recordKey) => JSON.stringify([recordKey, 'f']);
  const row = (recordKey, fields) => ({ recordKey, biologicalFeatureId: 'f', featureVisibility: null,
    labelVisibility: null, labelText: null, labelSourceText: null, ...fields });
  const edits = {
    featureOverrides: {
      [key('record-2')]: row('record-2', { labelText: 'own' }),
      [key('circular-B-#2')]: row('circular-B-#2', { labelText: 'alias', featureVisibility: 'off' })
    },
    featurePlacementOverrides: {}, featureStrokeOverrides: {}, featureColorOverrides: {}, annotationSets: []
  };
  assert.deepEqual(movedRecordEdits(edits, recordKeyMoves(canvas, canvas)).featureOverrides,
    { [key('record-2')]: row('record-2', { labelText: 'own', featureVisibility: 'off' }) });
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
