// @ts-check
// The twin registry (boundary audit M2): tools/twins.json lists every piece of
// logic that Python and the Web each implement, and `twin.unregistered.v1`
// (tools/web-twin-detectors.mjs) finds the twins it can see in source text.
// The registry and the detector agree exactly: a new twin needs a new row, and
// a row whose twin is gone leaves the registry with it.
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
  assert.deepEqual(check.unregistered, [], 'a twin without a row in tools/twins.json: delete one side, or add the row in an authority-only pull request');
  assert.deepEqual(check.stale, [], 'a row whose twin the detector no longer finds: remove the row');
});

test('every owner and vector a twin row names exists', () => {
  const rows = /** @type {{ twins: { id: string, pythonOwner: string[], jsOwner: string[], vector: string | null }[] }} */ (REGISTRY).twins;
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
  const [first] = /** @type {{ twins: { id: string }[] }} */ (REGISTRY).twins.filter(({ id }) => /^(?:vector|function):/.test(id));
  assert.ok(first, 'the registry has a detected row');
  const gone = checkTwinRegistry(REGISTRY, detectTwins(SOURCES).filter(({ id }) => id !== first.id));
  assert.deepEqual(gone.stale, [first.id]);
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
  assert.deepEqual(check.errors, [
    'tools/twins.json: twins[0].vector: expected tests/fixtures/<name>.json or null',
    'tools/twins.json: twins[1].id: duplicate vector:x',
    'tools/twins.json: twins[1].concern: expected text',
    'tools/twins.json: twins[1].pythonOwner: expected a list of paths',
    'tools/twins.json: twins[2].id: expected vector:, function:, constant: or concern: and a name',
    'tools/twins.json: twins[2].retireWhen: expected text'
  ]);
  assert.deepEqual(checkTwinRegistry([], []).errors, [
    'tools/twins.json: schemaVersion must be 1',
    'tools/twins.json: twins must be an array'
  ]);
});
