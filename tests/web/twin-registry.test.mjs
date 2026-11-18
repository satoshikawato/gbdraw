// @ts-check
// The twin registry (boundary audit M2): tools/twins.json lists every piece of
// logic that Python and the Web each implement, and `twin.unregistered.v1`
// (tools/web-twin-detectors.mjs) finds the twins it can see in source text.
// The registry and the detector agree exactly: a new twin needs a new row, and
// a row whose twin is gone leaves the registry with it. A new row goes in
// before its twin as pending, and the change that adds the twin removes it.
// A row the base lacks is admitted with its fields and paid for by retired
// rows or a waiver (`checkTwinAdditions`, run by the Gate against the base).
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import {
  checkTwinAdditions,
  checkTwinRegistry,
  detectTwins,
  readTwinExchangeBalance,
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
  // An owner `path:name` names a function its file still holds: a move of the
  // function leaves the row naming the old module (OV-2505).
  const unnamed = rows.flatMap(({ id, pythonOwner, jsOwner }) => [...pythonOwner, ...jsOwner]
    .flatMap((owner) => {
      const [, path, name] = /^(.*?)(?::([^/]*))?$/.exec(owner) ?? [];
      if (!name || !existsSync(join(REPOSITORY_ROOT, path))) return [];
      const word = name.split('.').pop()?.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') ?? '';
      return new RegExp(`\\b${word}\\b`).test(readFileSync(join(REPOSITORY_ROOT, path), 'utf8')) ? [] : [`${id}: ${owner}`];
    }));
  assert.deepEqual(unnamed, []);
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

test('the admission fields of a row are well formed wherever they appear', () => {
  const row = { id: 'concern:x', concern: 'c', pythonOwner: ['a.py'], jsOwner: ['a.js'], vector: null, reason: 'r', retireWhen: 'w', source: 's' };
  const registry = (extra) => ({ schemaVersion: 1, twins: [{ ...row, ...extra }] });
  assert.deepEqual(checkTwinRegistry(registry({ whyNotAskPython: 'p', whyNotGenerate: 'g', retireBy: '2027-01-31', waiver: 'w' }), []).errors, []);
  assert.deepEqual(checkTwinRegistry(registry({ whyNotAskPython: ' ', whyNotGenerate: 1, retireBy: '2027-02-30', waiver: '' }), []).errors, [
    'tools/twins.json: twins[0].whyNotAskPython: expected text',
    'tools/twins.json: twins[0].whyNotGenerate: expected text',
    'tools/twins.json: twins[0].waiver: expected text',
    'tools/twins.json: twins[0].retireBy: expected a date YYYY-MM-DD'
  ]);
  assert.deepEqual(checkTwinRegistry(registry({ retireBy: '2027-1-31' }), []).errors, [
    'tools/twins.json: twins[0].retireBy: expected a date YYYY-MM-DD'
  ]);
});

test('a row the bases lack needs its admission fields, and retired rows or a waiver pay for it', () => {
  const row = (id, extra = {}) => ({ id, concern: 'c', pythonOwner: ['a.py'], jsOwner: ['a.js'], vector: null, reason: 'r', retireWhen: 'w', source: 's', ...extra });
  const admitted = { whyNotAskPython: 'p', whyNotGenerate: 'g', retireBy: '2027-03-31' };
  const base = { schemaVersion: 1, twins: [row('concern:old')] };
  const head = (...rows) => ({ schemaVersion: 1, twins: [row('concern:old', { reason: 'edited' }), ...rows] });
  let reads = 0;
  // The day of the check is the UTC day: 2026-11-20 23:30 UTC is 2026-11-21 in Japan.
  const options = (balance) => ({ now: new Date('2026-11-20T23:30:00Z'), balance: () => { reads += 1; return balance; } });
  // A row the base has needs none of the fields, and no addition reads no balance.
  assert.deepEqual(checkTwinAdditions([base], head(), options(0)), []);
  assert.equal(reads, 0);
  assert.deepEqual(checkTwinAdditions([base], head(row('concern:new')), options(1)), [
    'tools/twins.json: concern:new: an added row needs whyNotAskPython, whyNotGenerate, retireBy (a date after 2026-11-20, UTC)'
  ]);
  assert.deepEqual(checkTwinAdditions([base], head(row('concern:new', { ...admitted, retireBy: '2026-11-20' })), options(1)), [
    'tools/twins.json: concern:new: an added row needs retireBy (a date after 2026-11-20, UTC)'
  ]);
  reads = 0;
  assert.deepEqual(checkTwinAdditions([base], head(row('concern:new', { ...admitted, retireBy: '2026-11-21' })), options(1)), []);
  assert.equal(reads, 1);
  // A row is added only when no base has it: the merge base may lack a row the trusted base has
  // (merged after the fork), and the trusted base may lack a row the merge base has (retired after it).
  reads = 0;
  const trusted = { schemaVersion: 1, twins: [row('concern:merged')] };
  assert.deepEqual(checkTwinAdditions([base, trusted], head(row('concern:merged')), options(0)), []);
  assert.deepEqual(checkTwinAdditions([trusted, base], head(row('concern:merged')), options(0)), []);
  assert.equal(reads, 0);
  const two = checkTwinAdditions([base], head(row('concern:new', admitted), row('concern:more', admitted)), options(1));
  assert.equal(two.length, 1);
  assert.match(two[0], /^tools\/twins\.json: 2 added row\(s\) without a waiver, but .* leave 1 to pay for additions: retire rows first, or give each added row a "waiver"/);
  // A row with a waiver takes nothing, and the balance is not read when every added row has one.
  reads = 0;
  assert.deepEqual(checkTwinAdditions([base], head(row('concern:new', { ...admitted, waiver: 'w' })), options(-3)), []);
  assert.equal(reads, 0);
  assert.deepEqual(checkTwinAdditions([base], head(row('concern:new', admitted), row('concern:waived', { ...admitted, waiver: 'w' })), options(1)), []);
  // Without a registry in the bases, every row is added.
  assert.equal(checkTwinAdditions([null, null], base, options(1)).length, 1);
});

test('the balance counts the rows production changes retired since the rule began, less the rows added without a waiver', () => {
  const root = mkdtempSync(join(tmpdir(), 'twin-balance-'));
  try {
    const environment = {
      ...process.env, GIT_CONFIG_GLOBAL: '/dev/null', GIT_CONFIG_NOSYSTEM: '1',
      GIT_AUTHOR_NAME: 't', GIT_AUTHOR_EMAIL: 't@example.com', GIT_COMMITTER_NAME: 't', GIT_COMMITTER_EMAIL: 't@example.com'
    };
    /** @param {string[]} args @param {string} date */
    const git = (args, date = '2026-11-20T10:00:00+09:00') => execFileSync('git', ['-C', root, ...args], {
      encoding: 'utf8', env: { ...environment, GIT_AUTHOR_DATE: date, GIT_COMMITTER_DATE: date }
    });
    /** @param {Record<string, string>} files */
    const write = (files) => Object.entries(files).forEach(([path, text]) => {
      mkdirSync(dirname(join(root, path)), { recursive: true });
      writeFileSync(join(root, path), text);
      git(['add', path]);
    });
    /** @param {string[]} ids @param {string[]} [waived] */
    const registry = (ids, waived = []) => JSON.stringify({
      schemaVersion: 1, twins: [...ids.map((id) => ({ id })), ...waived.map((id) => ({ id, waiver: 'w' }))]
    });
    /** @param {string} branch @param {Record<string, string>[]} commits @param {string} date */
    const merge = (branch, commits, date) => {
      git(['switch', '-q', '-c', branch]);
      commits.forEach((files) => {
        write(files);
        git(['commit', '-q', '-m', branch], date);
      });
      git(['switch', '-q', 'dev']);
      git(['merge', '-q', '--no-ff', branch, '-m', `Merge ${branch}`], date);
    };
    git(['init', '-q', '-b', 'dev']);
    write({ [TWIN_REGISTRY_PATH]: registry(['concern:a', 'concern:b', 'concern:c', 'concern:d', 'concern:e']), 'gbdraw/x.py': 'x = 0\n' });
    git(['commit', '-q', '-m', 'root'], '2026-11-10T10:00:00+09:00');
    // Before the rule began, a production change retires a row: it pays for nothing.
    merge('retire-e', [{ [TWIN_REGISTRY_PATH]: registry(['concern:a', 'concern:b', 'concern:c', 'concern:d']), 'gbdraw/x.py': 'x = 1\n' }], '2026-11-11T10:00:00+09:00');
    // A production change retires two rows, counted once through the first-parent merge.
    merge('retire-ab', [{ [TWIN_REGISTRY_PATH]: registry(['concern:c', 'concern:d']), 'gbdraw/web/js/services/a.js': 'export const a = 1;\n' }], '2026-11-18T10:00:00+09:00');
    // A change without production code removes a row: it pays for nothing.
    merge('drop-c', [{ [TWIN_REGISTRY_PATH]: registry(['concern:d']), 'tests/test_x.py': 'x = 1\n' }], '2026-11-19T10:00:00+09:00');
    // A merge dated before the rule began hides none of the changes behind it.
    merge('misdated', [{ 'docs/x.md': 'x\n' }], '2026-11-12T10:00:00+09:00');
    // An addition with a waiver takes nothing; one without takes one. The rows are read as the merge
    // brought them: the branch's first commit, before concern:g had its waiver, is not counted.
    merge('add-fg', [
      { [TWIN_REGISTRY_PATH]: registry(['concern:d', 'concern:f', 'concern:g']) },
      { [TWIN_REGISTRY_PATH]: registry(['concern:d', 'concern:f'], ['concern:g']) }
    ], '2026-11-20T10:00:00+09:00');
    /** @type {string[]} the git commands each reader ran */
    const commands = [];
    /** @param {string} repository @returns {(args: string[], input?: string) => string} */
    const reader = (repository) => (args, input) => {
      commands.push(args.find((arg, index) => !arg.startsWith('-') && args[index - 1] !== '-c') ?? '');
      return execFileSync('git', ['-C', repository, ...args], { encoding: 'utf8', env: environment, input, stdio: ['pipe', 'pipe', 'ignore'] });
    };
    /** @param {string} path */
    const production = (path) => path.startsWith('gbdraw/');
    assert.deepEqual(
      ['dev~4', 'dev~3', 'dev~2', 'dev~1', 'dev'].map((revision) => readTwinExchangeBalance(reader(root), revision, production)),
      [0, 2, 2, 2, 1]
    );
    // A clone without blobs (the trusted-base Gate) fetches the registry versions it reads; when it
    // cannot, the balance fails instead of reading the registry as absent.
    git(['config', 'uploadpack.allowFilter', 'true']);
    const clone = join(root, '.partial');
    execFileSync('git', ['clone', '-q', '--bare', '--filter=blob:none', `file://${root}`, clone], { env: environment, stdio: 'ignore' });
    commands.length = 0;
    assert.equal(readTwinExchangeBalance(reader(clone), 'dev', production), 1);
    // One fetch brings every registry version the balance reads, before cat-file reads them.
    assert.deepEqual(commands.filter((command) => command !== 'diff'), ['log', 'config', 'fetch', 'cat-file']);
    const unreachable = join(root, '.unreachable');
    execFileSync('git', ['clone', '-q', '--bare', '--filter=blob:none', `file://${root}`, unreachable], { env: environment, stdio: 'ignore' });
    execFileSync('git', ['-C', unreachable, 'remote', 'set-url', 'origin', `file://${root}/missing`], { env: environment });
    assert.throws(() => readTwinExchangeBalance(reader(unreachable), 'dev', production));
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
});
