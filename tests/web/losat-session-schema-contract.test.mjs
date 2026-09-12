import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../../', import.meta.url));
const source = path => readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');
const numericConstant = (text, name) => {
  const matches = [...text.matchAll(new RegExp(`\\bconst\\s+${name}\\s*=\\s*(\\d+)\\s*;`, 'g'))];
  assert.equal(matches.length, 1, `Expected one declaration of ${name}`);
  return Number(matches[0][1]);
};

// Load the fallback's declarations without selecting an adapter or a browser.
// Compare independent schema namespaces; released input fixtures keep their
// original versions and are not current-writer expectations.
const python = JSON.parse(execFileSync('python', ['-c', `
import json, runpy
from gbdraw.session_io import CURRENT_SESSION_VERSION, PROTEIN_LOSAT_CACHE_SCHEMA, LOSAT_DERIVED_CACHE_SCHEMA
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
adapter = runpy.run_path('tests/run_losat_cache_browser_acceptance.py')
print(json.dumps({
    'writer': {'session': CURRENT_SESSION_VERSION, 'request': CANONICAL_REQUEST_SCHEMA,
               'proteinRaw': PROTEIN_LOSAT_CACHE_SCHEMA, 'proteinDerived': LOSAT_DERIVED_CACHE_SCHEMA},
    'adapter': {'session': adapter['CURRENT_SESSION_VERSION'], 'request': adapter['CURRENT_RENDER_REQUEST_SCHEMA'],
                'proteinRaw': adapter['CURRENT_PROTEIN_RAW_SCHEMA'], 'proteinDerived': adapter['CURRENT_PROTEIN_DERIVED_SCHEMA']}
}))
`], { cwd: root, encoding: 'utf8' }));

test('Web and Python current Session, request and protein-cache schemas agree', () => {
  const cache = source('gbdraw/web/js/app/losat-cache.js');
  assert.deepEqual({
    session: numericConstant(source('gbdraw/web/js/services/config.js'), 'SESSION_VERSION'),
    request: numericConstant(source('gbdraw/web/js/services/session-request.js'), 'CANONICAL_REQUEST_SCHEMA'),
    proteinRaw: numericConstant(cache, 'PROTEIN_LOSAT_CACHE_SCHEMA'),
    proteinDerived: numericConstant(cache, 'LOSAT_DERIVED_CACHE_SCHEMA')
  }, python.writer);
});

test('Node LOSAT acceptance expects the current writer for every output schema', () => {
  const adapter = source('tests/web/losat-cache-migration.playwright.spec.js');
  assert.deepEqual({
    session: numericConstant(adapter, 'CURRENT_SESSION_VERSION'),
    request: numericConstant(adapter, 'CURRENT_RENDER_REQUEST_SCHEMA'),
    proteinRaw: numericConstant(adapter, 'CURRENT_PROTEIN_RAW_SCHEMA'),
    proteinDerived: numericConstant(adapter, 'CURRENT_PROTEIN_DERIVED_SCHEMA')
  }, python.writer);
});

test('Python LOSAT acceptance expects the current writer for every output schema', () => {
  assert.deepEqual(python.adapter, python.writer);
});
