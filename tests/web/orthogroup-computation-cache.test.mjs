import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
globalThis.window = { Vue: context.Vue, setTimeout, clearTimeout };
const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
const { rekeyOrthogroupOverrides } = await import('../../gbdraw/web/js/services/orthogroup-feature-metadata.js');
const { ref, reactive } = context.Vue;

test('group lookup reuses valid identities and invalidates on ID edits, new aliases, duplicates and replacement', () => {
  let reads = 0;
  const source = new Proxy({ id: 'one', members: [] }, {
    get(target, key, receiver) {
      if (key === 'id') reads += 1;
      return Reflect.get(target, key, receiver);
    }
  });
  const state = {
    orthogroups: ref([source]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('one'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]), biologicalFeatures: ref([])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const group = state.orthogroups.value[0];
  assert.equal(editor.getOrthogroupById('one'), group);
  const rows = editor.orthogroupRows.value;
  assert.equal(rows[0].name, 'one');
  const afterFirstLookup = reads;
  for (let i = 0; i < 20; i += 1) assert.equal(editor.getOrthogroupById('one'), group);
  state.showRightDrawer.value = true;
  assert.equal(editor.getOrthogroupById('one'), group);
  assert.equal(reads, afterFirstLookup, 'unrelated UI updates do not rescan group identities');
  assert.equal(editor.orthogroupRows.value, rows, 'unchanged row metadata is reused');
  state.orthogroupDescriptionOverrides.one = 'Edited description';
  assert.equal(editor.orthogroupRows.value[0].renamed, true);
  delete state.orthogroupDescriptionOverrides.one;
  assert.equal(editor.orthogroupRows.value[0].renamed, false);
  // CO-08 (D-21): a name belongs to the member set, not to the ID string. An
  // ID edit re-indexes the lookup; the name follows only by member set.
  state.orthogroupNameOverrides.one = 'User name';
  const before = [{ id: 'one', members: [{ proteinId: 'h_a' }] }];
  group.id = 'renamed';
  assert.equal(editor.getOrthogroupById('one'), null);
  assert.equal(editor.getOrthogroupById('renamed'), group);
  assert.equal(editor.orthogroupRows.value[0].name, 'renamed', 'an ID string alone carries no name');
  const followed = rekeyOrthogroupOverrides({
    previousGroups: before, candidateGroups: [{ id: 'renamed', members: [{ proteinId: 'h_a' }] }],
    names: state.orthogroupNameOverrides, descriptions: {}, dormant: {}
  });
  assert.deepEqual(followed, { names: { renamed: 'User name' }, descriptions: {}, dormant: {} });
  delete state.orthogroupNameOverrides.one;
  group.orthogroup_id = 'conflict';
  assert.equal(editor.getOrthogroupById('renamed'), null);
  delete group.orthogroup_id;
  assert.equal(editor.getOrthogroupById('renamed'), group);
  state.orthogroups.value.push({ orthogroupId: 'renamed', members: [] });
  assert.equal(editor.getOrthogroupById('renamed'), null);
  state.orthogroups.value.pop();
  assert.equal(editor.getOrthogroupById('renamed'), group);
  state.orthogroups.value = [{ id: 'replacement', members: [] }];
  assert.equal(editor.getOrthogroupById('renamed'), null);
  assert.equal(editor.getOrthogroupById('replacement'), state.orthogroups.value[0]);
});

test('group names follow exact member sets and unmatched names stay dormant (CO-08, D-21)', () => {
  const group = (id, ...members) => ({ id, members: members.map((proteinId) => ({ proteinId })) });
  const first = rekeyOrthogroupOverrides({
    previousGroups: [group('og_16', 'a', 'b'), group('og_17', 'c', 'd')],
    candidateGroups: [group('og_16', 'c', 'd', 'e'), group('og_17', 'a', 'b')],
    names: { og_16: 'LivB', og_17: 'LivC' },
    descriptions: { og_17: 'transporter' },
    dormant: {}
  });
  assert.deepEqual(first.names, { og_17: 'LivB' }, 'the name moves with {a,b} and never stays on og_16');
  assert.deepEqual(first.descriptions, {});
  assert.deepEqual(Object.values(first.dormant), [{ name: 'LivC', description: 'transporter' }]);
  const restored = rekeyOrthogroupOverrides({
    previousGroups: [group('og_16', 'c', 'd', 'e'), group('og_17', 'a', 'b')],
    candidateGroups: [group('og_3', 'c', 'd'), group('og_4', 'b', 'a')],
    names: first.names, descriptions: first.descriptions, dormant: first.dormant
  });
  assert.deepEqual(restored, {
    names: { og_3: 'LivC', og_4: 'LivB' }, descriptions: { og_3: 'transporter' }, dormant: {}
  }, 'a dormant name returns when its member set forms one group again');
});

test('member tables name each record in plain text, as the diagram heading reads (OV-369)', () => {
  const state = {
    orthogroups: ref([]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref(''), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), extractedFeatures: ref([]), biologicalFeatures: ref([]),
    linearSeqs: reactive([
      { uid: 'record-1', definition: '<i>Streptomyces lividus</i> CBS 844.73', gb: { name: 'BGC0000708' } },
      { uid: 'record-2', definition: '', file_definition: '<i>S. fradiae</i>', gb: { name: 'BGC0000709' } },
      { uid: 'record-3', definition: '', gb: { name: 'BGC0000711' } }
    ])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const rows = editor.groupOrthogroupMembersByRecord([
    { recordIndex: 0 }, { recordIndex: 1 }, { recordIndex: 2 }
  ]);
  assert.deepEqual(rows.map(({ recordLabel }) => recordLabel),
    ['Streptomyces lividus CBS 844.73', 'S. fradiae', 'BGC0000711']);
});
