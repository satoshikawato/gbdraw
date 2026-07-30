import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-definition-groups-'));
const files = [
  'app/results.js',
  'app/record-groups.js',
  'app/legend/utils.js',
  'app/legend-layout/transform-utils.js',
  'services/svg-serialization.js'
];
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');
await Promise.all(
  files.map(async (relativePath) => {
    const destination = join(tempDir, relativePath);
    await mkdir(dirname(destination), { recursive: true });
    await writeFile(
      destination,
      await readFile(join('gbdraw/web/js', relativePath), 'utf8'),
      'utf8'
    );
  })
);

const {
  findRecordDefinitionGroup,
  findSingleRecordDefinitionGroup,
  preserveDefinitionGroupDomIdentity
} = await import(pathToFileURL(join(tempDir, 'app/results.js')));

const group = ({
  id,
  recordId = null,
  recordIndex = null,
  role = 'record-definition',
  transform = null
}) => {
  const attributes = new Map([['id', id]]);
  if (recordId !== null) attributes.set('data-gbdraw-record-id', String(recordId));
  if (recordIndex !== null) {
    attributes.set('data-gbdraw-record-index', String(recordIndex));
  }
  if (role !== null) attributes.set('data-gbdraw-role', role);
  if (transform !== null) attributes.set('transform', transform);
  return {
    tagName: 'g',
    id,
    parentElement: null,
    getAttribute: (name) => attributes.get(name) ?? null,
    hasAttribute: (name) => attributes.has(name),
    setAttribute(name, value) {
      attributes.set(name, String(value));
      if (name === 'id') this.id = String(value);
    }
  };
};

const svg = (groups) => {
  const element = {
    tagName: 'svg',
    children: groups,
    parentElement: null,
    querySelectorAll(selector) {
      if (selector === 'g[data-gbdraw-role="record-definition"]') {
        return groups.filter(
          (candidate) =>
            candidate.getAttribute('data-gbdraw-role') === 'record-definition'
        );
      }
      if (selector === 'g[id$="_definition"]') {
        return groups.filter((candidate) => candidate.id.endsWith('_definition'));
      }
      return [];
    },
    getElementById(id) {
      return groups.find((candidate) => candidate.id === id) || null;
    }
  };
  groups.forEach((candidate) => {
    candidate.parentElement = element;
  });
  return element;
};

const firstDuplicate = group({
  id: 'same_definition_record_1',
  recordId: 'same',
  recordIndex: 0
});
const secondDuplicate = group({
  id: 'same_definition_record_2__record_2',
  recordId: 'same',
  recordIndex: 1,
  transform: 'translate(25, 40)'
});
const duplicateSvg = svg([firstDuplicate, secondDuplicate]);

assert.equal(
  findRecordDefinitionGroup(duplicateSvg, {
    definition_group_id: 'same_definition',
    record_index: 0
  }),
  firstDuplicate
);
assert.equal(
  findRecordDefinitionGroup(duplicateSvg, {
    definition_group_id: 'same_definition',
    record_index: 1
  }),
  secondDuplicate
);

const imported = group({
  id: 'same_definition_record_2',
  recordId: 'same',
  recordIndex: 1
});
preserveDefinitionGroupDomIdentity(secondDuplicate, imported);
assert.equal(imported.id, 'same_definition_record_2__record_2');
assert.equal(imported.getAttribute('transform'), 'translate(25, 40)');

const semanticSingle = group({
  id: 'opaque-generated-id',
  recordId: '123/unsafe id',
  recordIndex: 0
});
const legacySingle = group({
  id: 'legacy_definition',
  role: null
});
const singleSvg = svg([legacySingle, semanticSingle]);
assert.equal(findSingleRecordDefinitionGroup(singleSvg), semanticSingle);

const fallbackSvg = svg([legacySingle]);
assert.equal(
  findRecordDefinitionGroup(fallbackSvg, {
    definition_group_id: 'legacy_definition',
    record_index: 0
  }),
  legacySingle
);

console.log('definition group tests passed');
