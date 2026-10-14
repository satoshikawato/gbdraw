// readSessionDocument (E2b-1): the band reader drops a content part it cannot
// read with a row whose path is the part's key, and R4 (`restoresNothing`)
// refuses a Session whose every held content part was dropped.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import {
  readSessionDocument, restoresNothing, SESSION_TOP_LEVEL_AUTHORITY
} from '../../gbdraw/web/js/services/session-authority.js';

const fixture = (name) => {
  const bytes = readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url));
  return JSON.parse(name.endsWith('.gz') ? gunzipSync(bytes) : bytes);
};
const V46 = 'legend-renamed-rows.v46.gbdraw-session.json.gz';
const V44 = 'lambda_basic_linear.v44-schema8.gbdraw-session.json.gz';
const rows = (document) => readSessionDocument(document).repairs.map(({ code, path }) => `${code}@${path}`);

test('a content part of another kind is dropped with one row at its key', () => {
  for (const [name, expected] of [
    [V46, {
      resources: 'session-part-dropped', results: 'session-part-dropped', editorState: 'session-part-dropped',
      modes: 'session-part-dropped', otherModeResult: 'session-part-dropped', orthogroupState: 'session-field-dropped',
      files: 'session-field-dropped', losatCache: 'cache-entry-dropped', losatDerivedCache: 'cache-entry-dropped',
      proteinIdentityManifest: 'cache-entry-dropped'
    }],
    [V44, {
      resources: 'session-part-dropped', results: 'session-part-dropped', editorState: 'session-part-dropped',
      features: 'session-field-dropped', orthogroupState: 'session-field-dropped', files: 'session-field-dropped',
      modes: 'session-field-dropped', otherModeResult: 'session-field-dropped', losatCache: 'cache-entry-dropped',
      losatDerivedCache: 'cache-entry-dropped', proteinIdentityManifest: 'cache-entry-dropped'
    }]
  ]) {
    const base = fixture(name);
    assert.deepEqual(rows(base), [], name);
    for (const [key, code] of Object.entries(expected)) {
      assert.notEqual(SESSION_TOP_LEVEL_AUTHORITY[key].part, null, key);
      assert.deepEqual(rows({ ...structuredClone(base), [key]: true }), [`${code}@${key}`], `${name} ${key}`);
    }
  }
});

test('R4: a Session whose every held content part was dropped restores nothing', () => {
  const stored = fixture(V46);
  const held = Object.keys(stored).filter((key) => SESSION_TOP_LEVEL_AUTHORITY[key]?.part
    && ![undefined, null].includes(stored[key]) && !(typeof stored[key] === 'object' && Object.keys(stored[key]).length === 0));
  assert.ok(held.includes('results') && held.includes('modes'), held.join());
  const dropped = (keys) => keys.map((path) => ({ code: 'session-part-dropped', path, effect: 'dropped' }));
  assert.equal(restoresNothing(stored, dropped(held)), true);
  assert.equal(restoresNothing(stored, dropped(held.filter((key) => key !== 'modes'))), false);
  // A row inside a part does not drop the part.
  assert.equal(restoresNothing(stored, dropped(held.map((key) => `${key}[0]`))), false);
  assert.equal(restoresNothing(stored, []), false);
});
