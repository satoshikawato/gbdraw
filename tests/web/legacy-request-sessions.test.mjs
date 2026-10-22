// Request schemas 1 and 2 (Sessions 31-33): the stored requests of `main`
// history with their provenance, and the Load oracle that a new reading of
// them must keep (PR-D). The oracle is the dump of
// `helpers/load-oracle-probe.mjs` for each full fixture, recorded in
// `legacy-request-load.expected.json` as one SHA-256 per state area and saved
// part. The committed request is recorded apart: a new reading changes it.
// The record is the Load of PR-D: the one request reader lifting schemas 1
// and 2 and promoting them to schema 9, and Load's configuration check
// running the decoder's migration of a stored full config (OV-790). It equals
// the base Load with that migration added, but for the differences the
// hand-offs of PR-D steps D2 and D3a class (no rows for Python's "none" small
// tick interval, OV-791; schema-1 comparisons read as their schema-2
// re-saves; the committed request of schema 9).
//
//   UPDATE_LEGACY_LOAD_EXPECTATIONS=1 node --test tests/web/legacy-request-sessions.test.mjs
//
// rewrites the recorded file; review its diff. To see what changed in an
// area, run the probe on the fixture at both commits and diff the dumps.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync, writeFileSync } from 'node:fs';
import { execFile } from 'node:child_process';
import { promisify } from 'node:util';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';

import { readVectors } from './helpers/legacy-request-vectors.mjs';

const SESSIONS = new URL('../fixtures/sessions/', import.meta.url);
const EXPECTED = new URL('legacy-request-load.expected.json', SESSIONS);
const PROBE = fileURLToPath(new URL('helpers/load-oracle-probe.mjs', import.meta.url));
const sha256 = (value) => createHash('sha256').update(value).digest('hex');
const gunzipped = (name) => gunzipSync(readFileSync(new URL(name, SESSIONS)));
// The corpus of the shared vector.
const corpus = readVectors(new URL('../fixtures/legacy-render-request-vectors.json', import.meta.url));

// The full Session fixtures copied from `main` history (README.md): the
// Session, its request schema, and where the corpus holds its request.
const SCHEMA_1_FIXTURES = {
  'lambda_basic_linear.v31-schema1.gbdraw-session.json.gz': {
    version: 31, mode: 'linear', sha256: '15c8f32338453ede2c11e9976a7d972483483b1e2e85820d13398d3304326596',
    commit: '10d3a3d28b3c9faa42db01c8bd0bd36b9be8433c', path: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json'
  },
  'HmmtDNA_basic_circular.v32-schema1.gbdraw-session.json.gz': {
    version: 32, mode: 'circular', sha256: 'e565f29c358a0c372d967a29a2c58654c87417fa120800a90a05f6e1528b3429',
    commit: 'be13a5c030bbfba9e874cbaf22a74d08015532a9', path: 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'
  }
};
const FULL_FIXTURES = [
  ...Object.keys(SCHEMA_1_FIXTURES),
  ...corpus.fixtureReferences.map(({ path }) => path.split('/').pop())
];

test('the legacy request corpus holds each stored request once, with every place that stores it', () => {
  const { requests, fixtureReferences } = corpus;
  assert.deepEqual(
    [requests.length, requests.filter(({ schema }) => schema === 1).length, requests.flatMap(({ occurrences }) => occurrences).length],
    [23, 10, 47]
  );
  assert.equal(new Set(requests.map(({ id }) => id)).size, requests.length);
  assert.equal(new Set(requests.map(({ sha256: digest }) => digest)).size, requests.length);
  for (const { id, schema, mode, sha256: digest, occurrences, renderRequest } of requests) {
    assert.deepEqual([renderRequest.schema, renderRequest.mode], [schema, mode], id);
    assert.equal(sha256(JSON.stringify(renderRequest)), digest, id);
    assert.ok(occurrences.length > 0, id);
    for (const { commit, blob, path, sessionVersion } of occurrences) {
      assert.match(`${commit} ${blob}`, /^[0-9a-f]{40} [0-9a-f]{40}$/, id);
      assert.match(path, /^gbdraw\/web\/gallery\/sessions\/[^/]+\.gbdraw-session\.json(\.gz)?$/, id);
      assert.ok([31, 32, 33].includes(sessionVersion), id);
    }
  }
  // The committed (33, 2) fixtures, by reference.
  assert.equal(fixtureReferences.length, 4);
  for (const { path, sessionVersion, schema, mode, sha256: digest, corpusRequest } of fixtureReferences) {
    const session = JSON.parse(gunzipped(path.split('/').pop()));
    assert.deepEqual([session.version, session.renderRequest.schema, session.renderRequest.mode], [sessionVersion, schema, mode], path);
    assert.equal(sha256(JSON.stringify(session.renderRequest)), digest, path);
    assert.equal(corpusRequest, requests.find((request) => request.sha256 === digest)?.id ?? null, path);
  }
});

test('the schema-1 Session fixtures are the Gallery Sessions their provenance names', () => {
  for (const [name, { version, mode, sha256: digest, commit, path }] of Object.entries(SCHEMA_1_FIXTURES)) {
    const bytes = gunzipped(name);
    assert.equal(sha256(bytes), digest, name);
    const session = JSON.parse(bytes);
    assert.deepEqual([session.version, session.renderRequest.schema, session.renderRequest.mode], [version, 1, mode], name);
    const stored = corpus.requests.find((request) => request.sha256 === sha256(JSON.stringify(session.renderRequest)));
    assert.ok(stored?.occurrences.some((occurrence) => occurrence.commit === commit && occurrence.path === path), name);
  }
});

// The recorded form of one probe dump: what the Load and Save did, and one
// SHA-256 (16 hex digits) per state area and saved part; the requests apart.
const recorded = (dump) => ({
  oracle: {
    load: {
      status: dump.load.status,
      error: dump.load.error && { code: dump.load.error.code, operation: dump.load.error.operation },
      dialog: dump.load.dialogs.flatMap(({ items }) => items),
      notices: dump.load.notices
    },
    save: { status: dump.save.status, error: dump.save.error?.code ?? null },
    areas: Object.fromEntries(Object.entries(dump.areas).map(([area, digest]) => [area, digest.slice(0, 16)]))
  },
  requests: {
    committed: dump.committed?.request ?? null,
    saved: dump.save.request ?? null
  }
});
const runProbe = promisify(execFile);
const probe = async (name) => JSON.parse((await runProbe(
  process.execPath, [PROBE, fileURLToPath(new URL(name, SESSIONS))],
  { maxBuffer: 1 << 28, cwd: fileURLToPath(new URL('../../', import.meta.url)) }
)).stdout);

// The file keeps each distinct area map once (`areaSets`).
const packed = (records) => {
  const areaSets = {};
  const pack = ({ oracle, requests }) => {
    const key = sha256(JSON.stringify(oracle.areas)).slice(0, 12);
    areaSets[key] = oracle.areas;
    return { oracle: { ...oracle, areas: key }, requests };
  };
  const fixtures = Object.fromEntries(Object.entries(records).map(([name, record]) => [name, pack(record)]));
  return { fixtures, areaSets: Object.fromEntries(Object.entries(areaSets).sort()) };
};
const unpacked = ({ fixtures, areaSets }) => Object.fromEntries(Object.entries(fixtures).map(([name, record]) => (
  [name, { ...record, oracle: { ...record.oracle, areas: areaSets[record.oracle.areas] } }]
)));

test('Load and Save of each legacy full fixture are what the record holds', async () => {
  const actual = Object.fromEntries(await Promise.all(FULL_FIXTURES.map(async (name) => [name, recorded(await probe(name))])));
  if (process.env.UPDATE_LEGACY_LOAD_EXPECTATIONS === '1') {
    writeFileSync(EXPECTED, `${JSON.stringify(packed(actual), null, 1)}\n`);
    return;
  }
  const expected = unpacked(JSON.parse(readFileSync(EXPECTED, 'utf8')));
  assert.deepEqual(Object.keys(actual), Object.keys(expected));
  for (const name of FULL_FIXTURES) {
    const { oracle } = actual[name];
    const want = expected[name].oracle;
    const changed = Object.keys({ ...oracle.areas, ...want.areas }).filter((area) => oracle.areas[area] !== want.areas[area]);
    assert.deepEqual(changed, [], `${name}: areas changed`);
    assert.deepEqual(oracle, want, name);
  }
  // Recorded apart: the committed and saved requests (a new reading changes them).
  const requests = (records) => Object.fromEntries(Object.entries(records).map(([name, { requests: stored }]) => [name, stored]));
  assert.deepEqual(requests(actual), requests(expected));
});
