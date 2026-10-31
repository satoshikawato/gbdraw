import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../../', import.meta.url));
const source = path => readFileSync(new URL(`../../${path}`, import.meta.url), 'utf8');

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

test('Node LOSAT acceptance expects the current writer for every output schema', async () => {
  const adapter = source('tests/web/losat-cache-migration.playwright.spec.js');
  // The Session version and the request schema come from the generated
  // contract through the shared Playwright helper, not from a copy in the spec.
  for (const name of ['CURRENT_SESSION_VERSION', 'CURRENT_RENDER_REQUEST_SCHEMA']) {
    assert.doesNotMatch(adapter, new RegExp(`\\bconst\\s+${name}\\s*=`), `${name} is copied into the spec`);
  }
  assert.match(adapter, /CURRENT_REQUEST_SCHEMA: CURRENT_RENDER_REQUEST_SCHEMA,\s*CURRENT_SESSION_VERSION,[^}]*\} = require\('\.\/helpers\/app-lifecycle\.cjs'\)/);
  const helper = source('tests/web/helpers/app-lifecycle.cjs');
  assert.match(helper, /const CURRENT_REQUEST_SCHEMA = generatedConstant\('CANONICAL_REQUEST_SCHEMA'\);/);
  assert.match(helper, /const CURRENT_SESSION_VERSION = generatedConstant\('SESSION_VERSION'\);/);
  // The protein cache schemas come from the generated contract too.
  assert.match(adapter, /const CURRENT_PROTEIN_RAW_SCHEMA = generatedConstant\('PROTEIN_LOSAT_CACHE_SCHEMA'\);/);
  assert.match(adapter, /const CURRENT_PROTEIN_DERIVED_SCHEMA = generatedConstant\('LOSAT_DERIVED_CACHE_SCHEMA'\);/);
  const generated = await import('../../gbdraw/web/js/contracts/schema.generated.js');
  assert.deepEqual({
    session: generated.SESSION_VERSION,
    request: generated.CANONICAL_REQUEST_SCHEMA,
    proteinRaw: generated.PROTEIN_LOSAT_CACHE_SCHEMA,
    proteinDerived: generated.LOSAT_DERIVED_CACHE_SCHEMA
  }, python.writer);
});

test('Python LOSAT acceptance expects the current writer for every output schema', () => {
  assert.deepEqual(python.adapter, python.writer);
});
