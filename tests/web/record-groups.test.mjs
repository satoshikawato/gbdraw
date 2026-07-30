import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-record-groups-'));
const modulePath = join(tempDir, 'record-groups.mjs');
await writeFile(
  modulePath,
  await readFile('gbdraw/web/js/app/record-groups.js', 'utf8'),
  'utf8'
);
const {
  closestRecordGroup,
  directRecordGroups,
  isMultiRecordCanvasSvg
} = await import(pathToFileURL(modulePath));

const svg = (children = []) => {
  const element = { tagName: 'svg', id: '', children, parentElement: null };
  children.forEach((child) => {
    child.parentElement = element;
  });
  return element;
};
const group = (
  id,
  recordId = null,
  recordIndex = null,
  children = [],
  role = null
) => {
  const attributes = new Map();
  if (recordId !== null) attributes.set('data-gbdraw-record-id', String(recordId));
  if (recordIndex !== null) attributes.set('data-gbdraw-record-index', String(recordIndex));
  if (role !== null) attributes.set('data-gbdraw-role', String(role));
  const element = {
    tagName: 'g',
    id,
    children,
    parentElement: null,
    hasAttribute: (name) => attributes.has(name),
    getAttribute: (name) => attributes.get(name) ?? null
  };
  children.forEach((child) => {
    child.parentElement = element;
  });
  return element;
};
const leaf = (parentElement) => ({ tagName: 'path', parentElement });
const other = { tagName: 'defs', id: '', hasAttribute: () => false };

const single = svg([other, group('record_group_one', 'one', 0)]);
assert.equal(directRecordGroups(single).length, 1);
assert.equal(isMultiRecordCanvasSvg(single), false);

const multi = svg([
  group('record_group_one', 'one', 0),
  group('record_group_two', 'two', 1)
]);
assert.equal(directRecordGroups(multi).length, 2);
assert.equal(isMultiRecordCanvasSvg(multi), true);

const repeatedCircularFeatureSlots = svg([
  group('track_slot_inner', 'one', 0),
  group('track_slot_outer', 'one', 0)
]);
assert.equal(directRecordGroups(repeatedCircularFeatureSlots).length, 2);
assert.equal(isMultiRecordCanvasSvg(repeatedCircularFeatureSlots), false);

const duplicateRecordIdsAtDifferentIndexes = svg([
  group('record_group_first', 'same', 0),
  group('record_group_second', 'same', 1)
]);
assert.equal(isMultiRecordCanvasSvg(duplicateRecordIdsAtDifferentIndexes), true);

const semanticWithoutIndexes = svg([
  group('track_slot_inner', 'one'),
  group('track_slot_outer', 'one')
]);
assert.equal(isMultiRecordCanvasSvg(semanticWithoutIndexes), false);

const legacy = svg([group('record_one'), group('record_two')]);
assert.equal(directRecordGroups(legacy).length, 2);
assert.equal(isMultiRecordCanvasSvg(legacy), true);

const semantic = group('record_group_semantic', 'semantic', 0);
svg([semantic]);
assert.equal(closestRecordGroup(leaf(semantic)), semantic);

const nestedSemantic = group('record_group_nested', 'one', 0);
const circularOuter = group('record_0', 'one', 0, [nestedSemantic]);
const circularMultiSvg = svg([circularOuter, group('record_1', 'two', 1)]);
assert.equal(closestRecordGroup(leaf(nestedSemantic)), circularOuter);
assert.equal(isMultiRecordCanvasSvg(circularMultiSvg), true);

const persistedNestedSemantic = group('record_group_persisted_nested', 'one', 0);
const persistedCircularOuter = group(
  'record_0',
  null,
  null,
  [persistedNestedSemantic]
);
svg([persistedCircularOuter, group('record_1')]);
assert.equal(
  closestRecordGroup(leaf(persistedNestedSemantic)),
  persistedCircularOuter
);

const recordPrefixedDefinition = group('record_one_definition');
svg([recordPrefixedDefinition]);
assert.equal(closestRecordGroup(recordPrefixedDefinition), null);

const semanticDefinition = group(
  'record_generated_definition_record_2',
  'same',
  1,
  [],
  'record-definition'
);
svg([semanticDefinition]);
assert.equal(closestRecordGroup(semanticDefinition), null);

const legacyMultiDefinition = group('record_generated_definition_record_2');
svg([legacyMultiDefinition]);
assert.equal(closestRecordGroup(legacyMultiDefinition), null);

console.log('record group tests passed');
