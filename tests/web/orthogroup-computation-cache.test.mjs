import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
const { createOrthogroupEditor } = await import('../../gbdraw/web/js/app/orthogroups.js');
const { t } = await import('../../gbdraw/web/js/utils/i18n.js');
const { rekeyOrthogroupOverrides } = await import('../../gbdraw/web/js/services/orthogroup-feature-metadata.js');
const { ref, reactive } = Vue;

test('group lookup reuses valid identities and invalidates on ID edits, new aliases, duplicates and replacement', () => {
  let reads = 0;
  const source = new Proxy({ id: 'one', members: [] }, {
    get(target, key, receiver) {
      if (key === 'id') reads += 1;
      return Reflect.get(target, key, receiver);
    }
  });
  const state = {
    mode: 'linear', orthogroups: ref([source]),
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

// The rename tables are read by own key (OV-553). An own-key check alone is not
// a Vue dependency, so the first rename of a group must still reach the list.
test('the group list, a rendered name and the search follow the first rename of a group', () => {
  const state = {
    mode: 'linear', orthogroups: ref([{ id: 'og_1', name: 'LivB', description: 'Python text', members: [] }]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('og_1'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]), biologicalFeatures: ref([])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const rendered = Vue.computed(() => [
    editor.resolveOrthogroupName(editor.selectedOrthogroup.value),
    editor.resolveOrthogroupDescription(editor.selectedOrthogroup.value)
  ]);
  assert.deepEqual([editor.orthogroupRows.value[0].name, ...rendered.value], ['LivB', 'LivB', 'Python text']);
  // One write at a time: each of the two reads keeps its own dependency.
  editor.setOrthogroupDescriptionOverride('og_1', 'My text');
  assert.deepEqual(rendered.value, ['LivB', 'My text']);
  editor.setOrthogroupNameOverride('og_1', 'My transporter');
  assert.deepEqual(rendered.value, ['My transporter', 'My text']);
  assert.deepEqual([editor.orthogroupRows.value[0].name, editor.orthogroupRows.value[0].renamed], ['My transporter', true]);
  editor.resetOrthogroupRename('og_1');
  state.orthogroupSearch.value = 'transporter';
  assert.equal(editor.filteredOrthogroups.value.length, 0);
  editor.setOrthogroupNameOverride('og_1', 'My transporter');
  assert.equal(editor.filteredOrthogroups.value.length, 1, 'a search finds the group by its new name');
});

// OV-556: the "Edited" marker asks whether a rename table has the group's key.
// The marker must follow a first rename by itself, not through another read of its reader.
test('the Edited marker of a group follows its first name and its first description', () => {
  const state = {
    mode: 'linear', orthogroups: ref([{ id: 'og_1', name: 'LivB', members: [] }, { id: 'constructor', name: 'Own name', members: [] }]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('og_1'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]), biologicalFeatures: ref([])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const edited = Vue.computed(() => ['og_1', 'constructor'].map((id) => editor.isOrthogroupRenamed(id)));
  assert.deepEqual(edited.value, [false, false], 'a group id that names an inherited member is not edited');
  state.orthogroupNameOverrides.og_1 = 'My transporter';
  assert.deepEqual(edited.value, [true, false]);
  delete state.orthogroupNameOverrides.og_1;
  assert.deepEqual(edited.value, [false, false]);
  state.orthogroupDescriptionOverrides.og_1 = 'My text';
  assert.deepEqual(edited.value, [true, false]);
  state.orthogroupDescriptionOverrides.constructor = 'Own text';
  assert.deepEqual(edited.value, [true, true]);
  delete state.orthogroupDescriptionOverrides.og_1;
  assert.deepEqual(edited.value, [false, true]);
  // The list rows carry the same marker.
  assert.deepEqual(editor.orthogroupRows.value.map((row) => [row.id, row.renamed]), [['og_1', false], ['constructor', true]]);
  state.orthogroupDescriptionOverrides.og_1 = 'My text';
  const rows = editor.orthogroupRows.value;
  assert.deepEqual(rows.map((row) => row.renamed), [true, true]);
  // The marker asks for the key, not for its text: typing a description rebuilds no row.
  state.orthogroupDescriptionOverrides.og_1 = 'My longer text';
  assert.equal(editor.orthogroupRows.value, rows);
  assert.deepEqual(edited.value, [true, true]);
});

// OV-559: the Description field shows the user's own text of the selected group
// (not Python's description). A group id is data, so the field reads it by own
// key, through the editor. The field is the markup of index.html, compiled by Vue.
test('the Description field shows only the group\'s own edit and follows its first write (OV-559)', () => {
  const state = {
    mode: 'linear', orthogroups: ref([{ id: 'constructor', name: 'Own name', description: 'Python text', members: [] }]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('constructor'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]), biologicalFeatures: ref([])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const markup = html.match(/<textarea id="similarity-group-description"[^>]*><\/textarea>/)?.[0];
  assert.ok(markup, 'the Description field is found');
  const render = Vue.compile(markup);
  // What the composition root hands the template (app/app-setup.js).
  const bindings = {
    t,
    orthogroupDescriptionFieldText: editor.orthogroupDescriptionFieldText,
    setOrthogroupDescriptionOverride: editor.setOrthogroupDescriptionOverride,
    get selectedOrthogroup() { return editor.selectedOrthogroup.value; }
  };
  const field = Vue.computed(() => render(bindings, []).props.value);
  assert.equal(field.value, '', 'no edit: the field is empty, whatever the group id is');
  editor.setOrthogroupDescriptionOverride('constructor', 'Own text');
  assert.equal(field.value, 'Own text');
  editor.resetOrthogroupRename('constructor');
  assert.equal(field.value, '');
  // The reader, by row and by id.
  editor.setOrthogroupDescriptionOverride('constructor', 'Own text');
  assert.equal(editor.orthogroupDescriptionFieldText(editor.selectedOrthogroup.value), 'Own text');
  assert.equal(editor.orthogroupDescriptionFieldText('constructor'), 'Own text');
  assert.equal(editor.orthogroupDescriptionFieldText('toString'), '');
  assert.equal(editor.orthogroupDescriptionFieldText(null), '');
  // The template reads no rename table by key, and the root binds the reader.
  assert.doesNotMatch(html, /orthogroup(?:Name|Description)Overrides\s*\[/);
  const appSetup = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
  assert.match(appSetup, /orthogroupDescriptionFieldText: orthogroupActions\.orthogroupDescriptionFieldText,/);
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
    mode: 'linear', orthogroups: ref([]),
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
  // A member's record index names a record of the Result's request, whose
  // recordKey is the card uid (record selection).
  const getCommittedRequest = () => ({ records: state.linearSeqs.map(({ uid }) => ({ recordKey: uid })) });
  const editor = createOrthogroupEditor({ state: withDrawings(state), getCommittedRequest });
  const rows = editor.groupOrthogroupMembersByRecord([
    { recordIndex: 0 }, { recordIndex: 1 }, { recordIndex: 2 }
  ]);
  assert.deepEqual(rows.map(({ recordLabel }) => recordLabel),
    ['Streptomyces lividus CBS 844.73', 'S. fradiae', 'BGC0000711']);
});

// OV-920: the drawer shows the selected group's members, their sequences, and
// the group's copy menu. One detail per selected group is built when the
// group or a sequence source changes, not on every render that reads it.
test('the selected group detail is built once per change of its inputs (OV-920)', () => {
  let memberReads = 0;
  const counted = (member) => new Proxy(member, {
    get(target, key, receiver) {
      memberReads += 1;
      return Reflect.get(target, key, receiver);
    }
  });
  const state = {
    mode: 'linear', orthogroups: ref([{ id: 'og_1', members: [
      counted({ recordIndex: 0, stableFeatureSvgId: 'stable-a', start: 0, end: 30, strand: '+' }),
      counted({ recordIndex: 1, stableFeatureSvgId: 'stable-b', start: 3, end: 60, strand: '-' })
    ] }, { id: 'og_2', members: [] }]),
    orthogroupNameOverrides: reactive({}), orthogroupDescriptionOverrides: reactive({}),
    selectedOrthogroupId: ref('og_1'), orthogroupSearch: ref(''), orthogroupSortMode: ref('id'),
    clickedFeature: ref(null), showRightDrawer: ref(false), rightDrawerTab: ref('features'),
    svgContainer: ref(null), linearSeqs: reactive([]), extractedFeatures: ref([]),
    biologicalFeatures: ref([
      { record_idx: 0, stable_svg_id: 'stable-a', nucleotide_sequence: 'AAAA', amino_acid_sequence: 'MK' },
      { record_idx: 1, stable_svg_id: 'stable-b', nucleotide_sequence: 'CCCC' }
    ])
  };
  const editor = createOrthogroupEditor({ state: withDrawings(state) });
  const group = state.orthogroups.value[0];
  const detail = editor.selectedOrthogroupDetail.value;
  assert.equal(detail.group, group);
  assert.equal(detail.id, 'og_1');
  // The detail answers what the template asked the per-call helpers before.
  const expected = editor.groupOrthogroupMembersByRecord(editor.getEnrichedOrthogroupMembers(group));
  assert.deepEqual(detail.membersByRecord, expected);
  assert.deepEqual(detail.hasSequence, {
    nt: editor.getOrthogroupSequenceCount(group, 'nt') > 0,
    aa: editor.getOrthogroupSequenceCount(group, 'aa') > 0
  });
  const shown = detail.membersByRecord.flatMap(({ members }) => members);
  assert.deepEqual(shown.map((member) => [detail.hasMemberSequence(member, 'nt'), detail.hasMemberSequence(member, 'aa')]),
    shown.map((member) => [editor.hasOrthogroupMemberSequence(member, 'nt'), editor.hasOrthogroupMemberSequence(member, 'aa')]));
  assert.deepEqual(shown.map((member) => [detail.hasMemberSequence(member, 'nt'), detail.hasMemberSequence(member, 'aa')]),
    [[true, true], [true, false]]);

  const afterBuild = memberReads;
  for (let render = 0; render < 20; render += 1) {
    const again = editor.selectedOrthogroupDetail.value;
    assert.equal(again, detail);
    again.membersByRecord.forEach(({ members }) => members.forEach((member) => again.hasMemberSequence(member, 'nt')));
  }
  state.showRightDrawer.value = true;
  state.orthogroupSearch.value = 'og';
  assert.equal(editor.selectedOrthogroupDetail.value, detail, 'unrelated UI state reuses the detail');
  assert.equal(memberReads, afterBuild, 'no member is enriched again');

  // A new sequence source and another selected group build a new detail.
  state.biologicalFeatures.value = [{ record_idx: 0, stable_svg_id: 'stable-a', amino_acid_sequence: 'MK' }];
  assert.deepEqual(editor.selectedOrthogroupDetail.value.hasSequence, { nt: false, aa: true });
  state.selectedOrthogroupId.value = 'og_2';
  assert.equal(editor.selectedOrthogroupDetail.value.id, 'og_2');
  assert.deepEqual(editor.selectedOrthogroupDetail.value.membersByRecord, []);
  state.orthogroups.value = [];
  assert.equal(editor.selectedOrthogroupDetail.value, null);
});

// OV-920: the Similarity groups panel is part of the app template, which
// renders on every edit with the drawer closed too. Its bindings read the
// selected group's detail; no binding enriches members on each render.
test('the Similarity groups panel binds the selected detail, not enriching helpers (OV-920)', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf('id="right-drawer-panel-orthogroups"');
  const end = html.indexOf('data-similarity-group-actions', start);
  assert.ok(start > 0 && end > start, 'the panel markup is found');
  const panel = html.slice(start, end);
  for (const helper of ['hasOrthogroupSequence(', 'hasOrthogroupMemberSequence(', 'getEnrichedOrthogroupMembers(',
    'similarityAlignmentDrawerReferenceOptions(', 'similarityAlignmentDrawerDisabledReason(']) {
    assert.equal(panel.includes(helper), false, `${helper} is not called by the panel`);
  }
  assert.match(panel, /selectedOrthogroupDetail\.membersByRecord/);
  assert.match(panel, /similarityAlignmentDrawer\.options/);
});
