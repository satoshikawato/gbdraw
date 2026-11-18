// The shown name of a resource that stores none (OV-3660): the Web's `withoutResourceIdPrefix` and Python's
// `_without_resource_id_prefix` give each case of the shared vector its name (tools/twins.json,
// vector:tests/fixtures/resource-id-prefix-vectors.json; tests/test_input_summary.py reads it for Python).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';

import { withoutResourceIdPrefix } from '../../gbdraw/web/js/services/session-request.js';

/** @type {{ cases: { case: string, resourceId: string, name: string, shown: string }[] }} */
const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/resource-id-prefix-vectors.json', import.meta.url), 'utf8'));

test('the Web names a resource that stores no name as the shared vector says', () => {
  assert.ok(VECTORS.cases.length > 0);
  for (const { case: label, resourceId, name, shown } of VECTORS.cases) {
    assert.equal(withoutResourceIdPrefix(resourceId, name), shown, label);
  }
});
