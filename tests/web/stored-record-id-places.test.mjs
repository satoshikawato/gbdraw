// OV-3311: a stored Linear record whose selector is a record ID that a card's record text reads as a `#N`
// (`#2`, `#+2`) draws the record of that ID after Load. The card cannot hold that ID as text, so Load reads the
// record's place in its file from the record list (`readStoredRecordIdPlaces`, services/record-discovery.js) into
// the stored-leaf readings, and the projection gives the card that place (`#k`, the text the Record select gives
// such a record). A Session without such a record reads no file.
import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

const sourceRoot = join(process.cwd(), 'gbdraw', 'web', 'js');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-stored-record-id-places-'));
await cp(sourceRoot, join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
const moduleOf = (path) => import(pathToFileURL(join(tempRoot, 'js', ...path.split('/'))).href);
const { CANONICAL_REQUEST_SCHEMA, projectCanonicalSessionRequest } = await moduleOf('services/session-request.js');
const { readStoredRecordIdPlaces } = await moduleOf('services/record-discovery.js');
const { recordSelectorNumber } = await moduleOf('services/record-draw-selection.js');

// Records named `first`, `second`, `#2` and `#+2`, in this order: the ID `#2` is the third record and `#+2` the
// fourth, so reading either text as a `#N` draws another record.
const NAMES = ['first', 'second', '#2', '#+2'];
const genbank = NAMES.map((name, index) => [
  `LOCUS       ${name.padEnd(16)} ${10 + index} bp    DNA     linear   UNK 01-JAN-1980`,
  `DEFINITION  ${name}.`, `ACCESSION   ${name}`, `VERSION     ${name}`, 'FEATURES             Location/Qualifiers',
  'ORIGIN', `        1 ${'a'.repeat(10 + index)}`, '//', ''
].join('\n')).join('');
const resource = {
  kind: 'genbank', name: 'named.gb', type: 'text/plain', size: genbank.length, lastModified: 0,
  encoding: 'base64', data: Buffer.from(genbank).toString('base64')
};
const linearRecord = (recordKey, selector) => ({
  recordKey, cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'named' }, selector, windows: [],
  presentation: { label: null, subtitle: null, reverseComplement: false, gridRow: null, gridColumn: null },
  display: { isCircular: null, startCoordinate: null }
});
const request = (mode, records) => ({
  schema: CANONICAL_REQUEST_SCHEMA, mode, grouping: 'single', records,
  diagramOptions: { featurePlacements: [], featureOverrides: [] }, output: { prefix: 'ov3311' }
});
const countingFileOf = () => {
  const read = [];
  return {
    read,
    fileOf: (resourceId) => {
      read.push(resourceId);
      return resourceId === 'named' ? new File([genbank], 'named.gb', { type: 'text/plain' }) : null;
    }
  };
};
// The record a card's text names, as Generate reads it (`selectorPayload`): a `#N` counts records, any other text
// is a record ID.
const drawnName = (text) => {
  const number = recordSelectorNumber(text);
  return number === null ? text : NAMES[number - 1];
};

test('a stored record ID that reads as #N draws the record of that ID after Load', async () => {
  const stored = request('linear', [
    linearRecord('record-1', { kind: 'recordId', value: '#2' }),
    linearRecord('record-2', { kind: 'recordId', value: '#+2' }),
    linearRecord('record-3', { kind: 'recordId', value: 'first' })
  ]);
  const { read, fileOf } = countingFileOf();
  const leafReadings = await readStoredRecordIdPlaces([stored], fileOf);
  assert.deepEqual(read, ['named', 'named']);
  const projected = projectCanonicalSessionRequest({
    renderRequest: stored, resources: { named: resource }, initializeCliInputs: true, leafReadings
  });
  const cards = projected.files.linearSeqs.map((seq) => seq.region_record_id);
  assert.deepEqual(cards.map(drawnName), ['#2', '#+2', 'first']);
  assert.deepEqual(cards, ['#3', '#4', 'first']);
});

test('Load reads no file when no stored Linear record ID reads as #N', async () => {
  const { read, fileOf } = countingFileOf();
  const places = await readStoredRecordIdPlaces([
    request('linear', [
      linearRecord('record-1', { kind: 'recordId', value: 'first' }),
      linearRecord('record-2', { kind: 'recordIndex', index: 1 }),
      linearRecord('record-3', null),
      // A `#` text that is no `#N` stays a record ID as it is.
      linearRecord('record-4', { kind: 'recordId', value: '#+ 2' })
    ]),
    // A Circular request names its records by identity, not by card text.
    request('circular', [linearRecord('record-1', { kind: 'recordId', value: '#2' })]),
    null
  ], fileOf);
  assert.deepEqual(read, []);
  assert.equal(places.size, 0);
});

test('an ID no record or several records hold has no place', async () => {
  const twice = `${genbank}${genbank}`;
  const places = await readStoredRecordIdPlaces([
    request('linear', [
      linearRecord('record-1', { kind: 'recordId', value: '#9' }),
      { ...linearRecord('record-2', { kind: 'recordId', value: '#2' }), source: { kind: 'genbank', resourceId: 'twice' } }
    ])
  ], (resourceId) => new File([resourceId === 'twice' ? twice : genbank], `${resourceId}.gb`));
  assert.equal(places.size, 0);
});
