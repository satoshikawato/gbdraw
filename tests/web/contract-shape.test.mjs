// @ts-check
// The Web checker of a gbdraw.contract_shapes description (utils/contract-shape.js)
// against the vectors it shares with its Python twin
// (tests/fixtures/contract_shape_vectors.json, tests/test_contract_shape_vectors.py),
// the generated request description, and the stored schema-9 requests.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readdirSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

import { CANONICAL_REQUEST_DESCRIPTION } from '../../gbdraw/web/js/contracts/request-shape.generated.js';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { SHAPE_VIOLATION_REASONS, readSpelling, shapeViolations } from '../../gbdraw/web/js/utils/contract-shape.js';

/** @import { ContractShape, LiteralShape } from '../../gbdraw/web/js/utils/contract-shape.js' */
/**
 * @typedef {{ path: string, reason: string }} Row
 * @typedef {{ id: string, shape: ContractShape, value: unknown, violations: Row[] }} Vector
 * @typedef {{ id: string, shape: ContractShape, text: string, python: Row[], web: Row[] }} JsonDifference
 * @typedef {{ id: string, shape: LiteralShape, value: unknown, read: string | null }} Spelling
 * @typedef {{ id: string, type: string, value: unknown, violations: Row[] }} RequestTypeVector
 * @typedef {{ types: Record<string, ContractShape>, reasons: string[], vectors: Vector[], jsonDifferences: JsonDifference[], spellings: Spelling[], requestTypes: RequestTypeVector[] }} Vectors
 */

const REPO = new URL('../../', import.meta.url);
/** @type {Vectors} */
const VECTORS = JSON.parse(readFileSync(new URL('tests/fixtures/contract_shape_vectors.json', REPO), 'utf8'));

test('the Web checker reports each shared vector as Python does', () => {
  for (const { id, shape, value, violations } of VECTORS.vectors) {
    assert.deepEqual(shapeViolations(value, shape, VECTORS.types, 'value'), violations, id);
  }
});

test('both checkers report the same closed set of reasons', () => {
  assert.deepEqual([...SHAPE_VIOLATION_REASONS], VECTORS.reasons);
  const reported = new Set(VECTORS.vectors.flatMap(({ violations }) => violations.map(({ reason }) => reason)));
  assert.deepEqual([...reported].sort(), [...SHAPE_VIOLATION_REASONS].sort());
});

test('the Web reads each listed JSON difference as listed', () => {
  for (const { id, shape, text, web } of VECTORS.jsonDifferences) {
    assert.deepEqual(shapeViolations(JSON.parse(text), shape, VECTORS.types, 'value'), web, id);
  }
});

test('the Web reads each shared spelling as Python does', () => {
  for (const { id, shape, value, read } of VECTORS.spellings) {
    assert.equal(readSpelling(value, shape), read, id);
  }
});

test('the generated request description reports each shared request-type vector as Python does', () => {
  for (const { id, type, value, violations } of VECTORS.requestTypes) {
    assert.deepEqual(
      shapeViolations(value, { kind: 'ref', name: type }, CANONICAL_REQUEST_DESCRIPTION.types, 'value'), violations, id
    );
  }
});

test('a violation holds a path and a reason, never the value', () => {
  const secret = 'private-genome-name';
  const found = shapeViolations({ a: secret, extra: secret }, {
    kind: 'object', fields: { a: { kind: 'integer' }, b: { kind: 'string' } }, optional: [], open: false
  }, {}, 'value');
  assert.deepEqual(found, [
    { path: 'value.extra', reason: 'unknown-key' },
    { path: 'value.b', reason: 'missing-key' },
    { path: 'value.a', reason: 'wrong-kind' }
  ]);
  assert.equal(JSON.stringify(found).includes(secret), false);
});

/**
 * Each shape of a description, with the name of the type it is in.
 * @param {Readonly<Record<string, ContractShape>>} types
 * @returns {[string, ContractShape][]}
 */
const allShapes = (types) => {
  /** @type {[string, ContractShape][]} */
  const found = [];
  /** @param {string} name @param {ContractShape} shape */
  const walk = (name, shape) => {
    found.push([name, shape]);
    if (shape.kind === 'array') walk(name, shape.items);
    else if (shape.kind === 'tuple') shape.items.forEach((item) => walk(name, item));
    else if (shape.kind === 'record') walk(name, shape.values);
    else if (shape.kind === 'union') shape.options.forEach((option) => walk(name, option));
    else if (shape.kind === 'object') Object.values(shape.fields).forEach((field) => walk(name, field));
  };
  Object.entries(types).forEach(([name, shape]) => walk(name, shape));
  return found;
};

test('the generated request description is the data Python describes', () => {
  const python = execFileSync('python', ['-c', [
    'import json',
    'from gbdraw.session_request_codec import describe_canonical_request',
    'print(json.dumps(describe_canonical_request()))'
  ].join('\n')], { cwd: fileURLToPath(REPO), encoding: 'utf8', maxBuffer: 1 << 26 });
  assert.deepEqual(CANONICAL_REQUEST_DESCRIPTION, JSON.parse(python));
});

test('the generated request description is frozen and names every shape it refers to', () => {
  const { root, types } = CANONICAL_REQUEST_DESCRIPTION;
  assert.ok(Object.hasOwn(types, root));
  /** @param {unknown} value @returns {boolean} */
  const deeplyFrozen = (value) => typeof value !== 'object' || value === null
    || (Object.isFrozen(value) && Object.values(value).every(deeplyFrozen));
  assert.ok(deeplyFrozen(CANONICAL_REQUEST_DESCRIPTION));
  const missing = allShapes(types).filter(([, shape]) => shape.kind === 'ref' && !Object.hasOwn(types, shape.name));
  assert.deepEqual(missing, []);
});

test('each default of the request description fits its field', () => {
  const { types } = CANONICAL_REQUEST_DESCRIPTION;
  /** @type {string[]} */
  const defaulted = [];
  for (const [name, shape] of allShapes(types)) {
    if (shape.kind !== 'object' || !shape.defaults) continue;
    for (const [field, value] of Object.entries(shape.defaults)) {
      defaulted.push(`${name}.${field}`);
      assert.deepEqual(shapeViolations(value, shape.fields[field], types, `${name}.${field}`), []);
    }
  }
  // The defaults a repair of a stored request needs at a part (OV-774).
  for (const field of [
    'CanonicalProteinComparisonSettings.collinearInferOrthogroups', 'CanonicalLinearLayout.recordGapPx',
    'CanonicalLinearLayout.similarityAlignment', 'CanonicalRecord.display',
    'CanonicalCircularRenderRequest.diagramOptions', 'CanonicalCircularRenderRequest.output',
    'CanonicalLinearRenderRequest.diagramOptions', 'CanonicalLinearRenderRequest.output'
  ]) {
    assert.ok(defaulted.includes(field), `${field} in ${defaulted.join(', ')}`);
  }
});

/**
 * The renderRequest of a Session 49's drawings, and of a 46-48 Session and of
 * its other mode's Result.
 * @param {unknown} session
 * @returns {[string, unknown][]}
 */
const storedRequests = (session) => {
  if (typeof session !== 'object' || session === null) return [];
  const { renderRequest, otherModeResult, drawings } = /** @type {{ renderRequest?: unknown, otherModeResult?: unknown, drawings?: { id: string, renderRequest?: unknown }[] }} */ (session);
  if (Array.isArray(drawings)) return drawings.map((drawing) => [`drawings[${drawing.id}].renderRequest`, drawing.renderRequest]);
  const other = typeof otherModeResult === 'object' && otherModeResult !== null
    ? /** @type {{ renderRequest?: unknown }} */ (otherModeResult).renderRequest
    : undefined;
  return /** @type {[string, unknown][]} */ ([['renderRequest', renderRequest], ['otherModeResult.renderRequest', other]]);
};

test('every stored schema-9 request of the Session fixtures and the Gallery fits the description', () => {
  const { root, types } = CANONICAL_REQUEST_DESCRIPTION;
  /** @type {Record<string, unknown>} */
  const problems = {};
  let checked = 0;
  for (const folder of ['tests/fixtures/sessions/', 'gbdraw/web/gallery/sessions/']) {
    for (const name of readdirSync(new URL(folder, REPO)).filter((file) => /\.json(?:\.gz)?$/.test(file)).sort()) {
      const raw = readFileSync(new URL(folder + name, REPO));
      const session = JSON.parse((name.endsWith('.gz') ? gunzipSync(raw) : raw).toString('utf8'));
      for (const [where, request] of storedRequests(session)) {
        if (/** @type {{ schema?: unknown } | null | undefined} */ (request)?.schema !== CANONICAL_REQUEST_SCHEMA) continue;
        checked += 1;
        const found = shapeViolations(request, types[root], types, where);
        if (found.length) problems[folder + name] = found;
      }
    }
  }
  assert.deepEqual(problems, {});
  assert.ok(checked >= 12, `${checked} schema-9 requests`);
});
