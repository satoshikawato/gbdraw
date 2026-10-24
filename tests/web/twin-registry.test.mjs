// @ts-check
// The twin registry (boundary audit M2): tools/twins.json lists every piece of
// logic that Python and the Web each implement, and `twin.unregistered.v1`
// (tools/web-twin-detectors.mjs) finds the twins it can see in source text.
// The registry and the detector agree exactly: a new twin needs a new row, and
// a row whose twin is gone leaves the registry with it. A new row goes in
// before its twin as pending, and the change that adds the twin removes it.
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import {
  checkTwinRegistry,
  detectTwins,
  readTwinSources,
  TWIN_DETECTOR_ID,
  TWIN_REGISTRY_PATH
} from '../../tools/web-twin-detectors.mjs';

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
/** @type {unknown} */
const REGISTRY = JSON.parse(readFileSync(join(REPOSITORY_ROOT, TWIN_REGISTRY_PATH), 'utf8'));
const SOURCES = readTwinSources(REPOSITORY_ROOT);

/**
 * @param {ReadonlyMap<string, string>} map
 * @param {Record<string, string>} added
 */
const withFiles = (map, added) => new Map([...map, ...Object.entries(added)]);

test(`${TWIN_DETECTOR_ID} reports exactly the registered twins`, () => {
  const check = checkTwinRegistry(REGISTRY, detectTwins(SOURCES));
  assert.deepEqual(check.errors, []);
  assert.deepEqual(check.unregistered, [], 'a twin without a row in tools/twins.json: delete one side, or add the row with "pending": true in an authority-only pull request first');
  assert.deepEqual(check.stale, [], 'a row whose twin the detector no longer finds: remove the row');
});

test('a vector in a subfolder of tests/fixtures is read, detected, and accepted as a row vector', () => {
  const vector = 'tests/fixtures/sessions/sub_twin_vectors.json';
  assert.ok(SOURCES.fixtures.some((path) => path.split('/').length > 3), 'readTwinSources lists subfolder vectors');
  const found = detectTwins({
    ...SOURCES,
    fixtures: [vector],
    tests: new Map([
      ['tests/web/sub-twin.test.mjs', "load('sub_twin_vectors.json');\n"],
      ['tests/test_sub_twin.py', 'load("sub_twin_vectors.json")\n']
    ])
  });
  const detected = found.filter(({ id }) => id.startsWith('vector:'));
  assert.deepEqual(detected.map(({ id }) => id), [`vector:${vector}`]);
  const row = { id: `vector:${vector}`, concern: 'c', pythonOwner: ['a.py'], jsOwner: ['a.js'], vector, reason: 'r', retireWhen: 'w', source: 's' };
  assert.deepEqual(checkTwinRegistry({ schemaVersion: 1, twins: [row] }, detected).errors, []);
});

test('every owner and vector a twin row that is not pending names exists', () => {
  const rows = /** @type {{ twins: { id: string, pythonOwner: string[], jsOwner: string[], vector: string | null, pending?: true }[] }} */ (REGISTRY).twins
    .filter(({ pending }) => pending !== true);
  const missing = rows.flatMap(({ id, pythonOwner, jsOwner, vector }) => [...pythonOwner, ...jsOwner, ...(vector ? [vector] : [])]
    .map((owner) => owner.replace(/:[^/]*$/, ''))
    .filter((path) => !existsSync(join(REPOSITORY_ROOT, path)))
    .map((path) => `${id}: ${path}`));
  assert.deepEqual(missing, []);
});

test('a twin the registry does not list fails, and so does a row without its twin', () => {
  const vector = 'tests/fixtures/extra_twin_vectors.json';
  const extra = detectTwins({
    ...SOURCES,
    fixtures: [...SOURCES.fixtures, vector],
    tests: withFiles(SOURCES.tests, {
      'tests/web/extra-twin.test.mjs': "const CASES = 'extra_twin_vectors.json';\n",
      'tests/test_extra_twin.py': 'CASES = FIXTURES / "extra_twin_vectors.json"\n'
    }),
    python: withFiles(SOURCES.python, { 'gbdraw/extra_twin.py': 'def migrate_extra_twin(value):\n    return value\n' }),
    services: withFiles(SOURCES.services, {
      'gbdraw/web/js/services/extra-twin.js': 'export const migrateExtraTwin = (value) => value;\n'
    })
  });
  const check = checkTwinRegistry(REGISTRY, extra);
  assert.deepEqual(check.unregistered, ['function:migrateExtraTwin', `vector:${vector}`]);
  assert.deepEqual(check.stale, []);
  const [first] = /** @type {{ twins: { id: string, pending?: true }[] }} */ (REGISTRY).twins
    .filter(({ id, pending }) => !pending && /^(?:vector|function):/.test(id));
  assert.ok(first, 'the registry has a detected row');
  const gone = checkTwinRegistry(REGISTRY, detectTwins(SOURCES).filter(({ id }) => id !== first.id));
  assert.deepEqual(gone.stale, [first.id]);
});

test('a pending row waits for its twin, and keeps pending only until the twin is detected', () => {
  const pendingRow = {
    id: 'function:migrateExtraTwin',
    concern: 'c',
    pythonOwner: ['gbdraw/extra_twin.py'],
    jsOwner: ['gbdraw/web/js/services/extra-twin.js'],
    vector: null,
    reason: 'r',
    retireWhen: 'w',
    source: 's',
    pending: true
  };
  const registry = { schemaVersion: 1, twins: [pendingRow] };
  assert.deepEqual(checkTwinRegistry(registry, []), { errors: [], unregistered: [], stale: [] });
  const twin = { id: pendingRow.id, python: ['gbdraw/extra_twin.py:migrate_extra_twin'], javascript: ['gbdraw/web/js/services/extra-twin.js:migrateExtraTwin'] };
  assert.deepEqual(checkTwinRegistry(registry, [twin]).errors, [
    'tools/twins.json: function:migrateExtraTwin: the twin is detected; remove pending'
  ]);
  const { pending: _pending, ...added } = pendingRow;
  assert.deepEqual(checkTwinRegistry({ schemaVersion: 1, twins: [added] }, [twin]), { errors: [], unregistered: [], stale: [] });
  assert.deepEqual(checkTwinRegistry({
    schemaVersion: 1,
    twins: [{ ...pendingRow, pending: false }, { ...pendingRow, id: 'constant:extra', pending: true }]
  }, []).errors, [
    'tools/twins.json: twins[0].pending: expected true on a vector: or function: row, or no pending',
    'tools/twins.json: twins[1].pending: expected true on a vector: or function: row, or no pending'
  ]);
});

test('the detector reads names, not comments, and compares Python and JavaScript spellings', () => {
  const detected = detectTwins({
    fixtures: ['tests/fixtures/a_cases.json', 'tests/fixtures/b_cases.json'],
    tests: new Map([
      ['tests/web/a.test.mjs', "// reads b_cases.json in Python only\nload('a_cases.json');\n"],
      ['tests/test_a.py', '# a_cases.json is read in JavaScript too\nload("a_cases.json")\nload("b_cases.json")\n'],
      ['tests/web/c.test.mjs', "load('xb_cases.json');\n"]
    ]),
    python: new Map([[
      'gbdraw/a.py',
      'def _split_draft_into_modes(draft):\n    """def normalize_hidden(x): in a docstring"""\n\ndef normalize(value):\n    pass\n\ndef resolve_label_color(value):\n    pass\n'
    ]]),
    services: new Map([[
      'gbdraw/web/js/services/a.js',
      [
        'export function splitDraftIntoModes(draft) { return draft; }',
        '// const normalizeHidden = (x) => x;',
        'const normalize = (value) => value;',
        'const resolveLabelColor = async (value) => value;',
        "const text = 'function readColorTable() {}';"
      ].join('\n')
    ]])
  });
  assert.deepEqual(detected, [
    { id: 'function:resolveLabelColor', python: ['gbdraw/a.py:resolve_label_color'], javascript: ['gbdraw/web/js/services/a.js:resolveLabelColor'] },
    { id: 'function:splitDraftIntoModes', python: ['gbdraw/a.py:_split_draft_into_modes'], javascript: ['gbdraw/web/js/services/a.js:splitDraftIntoModes'] },
    { id: 'vector:tests/fixtures/a_cases.json', python: ['tests/test_a.py'], javascript: ['tests/web/a.test.mjs'] }
  ]);
});

test('a malformed twin row is reported', () => {
  const check = checkTwinRegistry({
    schemaVersion: 1,
    twins: [
      { id: 'vector:x', concern: 'c', pythonOwner: [], jsOwner: [], vector: 'fixtures/x.json', reason: 'r', retireWhen: 'w', source: 's' },
      { id: 'vector:x', concern: '', pythonOwner: 'gbdraw/a.py', jsOwner: [], vector: null, reason: 'r', retireWhen: 'w', source: 's' },
      { id: 'twin:x', concern: 'c', pythonOwner: [], jsOwner: [], vector: null, reason: 'r', retireWhen: ' ', source: 's' }
    ]
  }, []);
  // A twin has two sides: a row with an empty owner list is not a twin.
  assert.deepEqual(check.errors, [
    'tools/twins.json: twins[0].pythonOwner: expected a non-empty list of paths',
    'tools/twins.json: twins[0].jsOwner: expected a non-empty list of paths',
    'tools/twins.json: twins[0].vector: expected tests/fixtures/<path>.json or null',
    'tools/twins.json: twins[1].id: duplicate vector:x',
    'tools/twins.json: twins[1].concern: expected text',
    'tools/twins.json: twins[1].pythonOwner: expected a non-empty list of paths',
    'tools/twins.json: twins[1].jsOwner: expected a non-empty list of paths',
    'tools/twins.json: twins[2].id: expected vector:, function:, constant: or concern: and a name',
    'tools/twins.json: twins[2].retireWhen: expected text',
    'tools/twins.json: twins[2].pythonOwner: expected a non-empty list of paths',
    'tools/twins.json: twins[2].jsOwner: expected a non-empty list of paths'
  ]);
  assert.deepEqual(checkTwinRegistry([], []).errors, [
    'tools/twins.json: schemaVersion must be 1',
    'tools/twins.json: twins must be an array'
  ]);
});
