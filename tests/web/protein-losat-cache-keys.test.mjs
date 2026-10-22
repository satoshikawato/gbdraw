import assert from 'node:assert/strict';
import test from 'node:test';

import { requireProteinLosatCacheKeys } from '../../gbdraw/web/js/app/run-analysis.js';

const refused = { name: 'Error', message: 'Protein cache key generation failed.' };

// The reply of the `buildProteinLosatCacheKeys` helper is unchecked data: Generate uses each key as a raw
// LOSAT cache key, so only one SHA-256 hex digest per pair is accepted.
test('a protein cache key reply holds one hex digest per pair, as text (OV-619)', () => {
  const key = 'a'.repeat(64);
  const other = 'b'.repeat(64);
  assert.deepEqual(requireProteinLosatCacheKeys({ keys: [key, other] }, 2), { keys: [key, other], derivedKey: '' });
  for (const keys of [[[key]], [new String(key)], [{ toString: () => key }]]) {
    assert.throws(() => requireProteinLosatCacheKeys({ keys }, 1), refused, JSON.stringify(keys));
  }
  for (const reply of [null, undefined]) {
    assert.throws(() => requireProteinLosatCacheKeys(reply, 1), refused, String(reply));
  }
});

test('a protein cache key reply of another shape is refused, and a reported error is thrown as it is', () => {
  const key = 'a'.repeat(64);
  for (const reply of [{}, { keys: null }, { keys: key }, { keys: [] }, { keys: [key, key] }, { keys: [key.slice(1)] },
    { keys: [key.toUpperCase()] }, { keys: [`${key}\n`] }, { keys: [5] }, 'keys', 5]) {
    assert.throws(() => requireProteinLosatCacheKeys(reply, 1), refused, JSON.stringify(reply));
  }
  const error = { code: 'RESOURCE_INVALID' };
  assert.throws(() => requireProteinLosatCacheKeys({ error, keys: [key] }, 1), (thrown) => thrown === error);
});

// OV-1063: the reply also holds the key of the derived payload when the request asked for it, and Generate
// looks the payload up by it, so it is one hex digest too.
test('a protein cache key reply holds the derived key the request asked for', () => {
  const key = 'a'.repeat(64);
  const derivedKey = 'c'.repeat(64);
  assert.deepEqual(
    requireProteinLosatCacheKeys({ keys: [key], derived: { key: derivedKey } }, 1, { derived: true }),
    { keys: [key], derivedKey }
  );
  assert.deepEqual(
    requireProteinLosatCacheKeys({ keys: [key], derived: { key: 'not read' } }, 1),
    { keys: [key], derivedKey: '' },
    'a request without derived inputs reads no derived key'
  );
  for (const derived of [undefined, null, {}, { key: null }, { key: key.slice(1) }, { key: key.toUpperCase() },
    { key: new String(derivedKey) }, { key: [derivedKey] }, derivedKey]) {
    assert.throws(
      () => requireProteinLosatCacheKeys({ keys: [key], derived }, 1, { derived: true }),
      refused,
      JSON.stringify(derived)
    );
  }
});
