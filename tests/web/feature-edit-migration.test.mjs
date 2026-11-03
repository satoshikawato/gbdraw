// Session 41-44 placement drafts: the JS reader that remains after Python took over the
// rendered-ID edit migration (tests/test_session_compat.py runs those vectors).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { migrateSessionFeaturePlacements } from '../../gbdraw/web/js/services/feature-edit-migration.js';

// The same vector pins the Python reader (tests/test_session_compat.py).
test('Session 41-44 placement drafts migrate to the vector both readers share', () => {
  const vector = JSON.parse(readFileSync(new URL('../fixtures/feature-placement-migration.json', import.meta.url), 'utf8'));
  assert.deepEqual(migrateSessionFeaturePlacements(structuredClone(vector.input)), vector.expected);
});
