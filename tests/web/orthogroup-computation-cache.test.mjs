import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import test from 'node:test';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
globalThis.window = { Vue: context.Vue, setTimeout, clearTimeout };
const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
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
  const editor = createOrthogroupEditor({ state });
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
  group.id = 'renamed';
  assert.equal(editor.getOrthogroupById('one'), null);
  assert.equal(editor.getOrthogroupById('renamed'), group);
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
