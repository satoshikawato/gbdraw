import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import { requireProteinLosatCacheKeys } from '../../gbdraw/web/js/app/run-analysis.js';
import { emptyProteinIdentityManifest, mergeProteinIdentityManifests } from '../../gbdraw/web/js/services/losat-cache.js';

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

// OV-1553: Generate asks the Worker for the keys of a pair list once per merged manifest and request. The block
// runs as written in run-analysis.js, with a Worker stand-in that counts the requests.
test('Generate asks for protein cache keys once per manifest and key request (OV-1553)', async () => {
  const runAnalysisPath = resolve('gbdraw/web/js/app/run-analysis.js');
  const source = await readFile(runAnalysisPath, 'utf8');
  const start = source.indexOf('          if (useProteinBlastp && preparedJobs.length > 0) {\n            throwIfGenerationCanceled();\n            const derivedRecords = [];');
  const end = source.indexOf('\n\n          // ', start);
  assert.ok(start > 0 && end > start);
  // The module's own cache and its readers, exported from a copy of the module.
  const probeDir = await mkdtemp(join(tmpdir(), 'gbdraw-protein-keys-'));
  const probePath = join(probeDir, 'run-analysis-keys-probe.mjs');
  const exported = ['getCachedFileRead', 'setCachedFileRead', 'proteinCacheKeysCache']
    .filter((name) => source.includes(`const ${name} =`));
  await writeFile(probePath, `${source.replace(/from '(\.\.?\/[^']+)'/g, (_match, specifier) => (
    `from '${new URL(specifier, pathToFileURL(runAnalysisPath)).href}'`
  ))}\nexport { ${exported.join(', ')} };\n`, 'utf8');
  const probe = await import(pathToFileURL(probePath));
  const requests = [];
  const block = new Function('context', `return (async () => {
    const {${['useProteinBlastp', 'preparedJobs', 'throwIfGenerationCanceled', 'useDerivedProteinPayloadCache',
      'converterRecordIndexes', 'getSeqEntry', 'getViewTransform', 'runDiagramHelperOperation',
      'DIAGRAM_HELPER_OPERATIONS', 'cloneJsonData', 'workingProteinIdentityManifest', 'normalizeLosatArgs',
      'derivedSettings', 'pairNumbers', 'requireProteinLosatCacheKeys', 'getCachedFileRead', 'setCachedFileRead',
      'proteinCacheKeysCache'].join(',')}} = context;
    let proteinCacheKeys = null;
    let derivedCacheKey = '';
    ${source.slice(start, end)}
    return { proteinCacheKeys, derivedCacheKey };
  })()`);
  const pair = (args) => ({ cacheMetadata: {
    queryRecordInstanceKey: 'r1', subjectRecordInstanceKey: 'r2', program: 'blastp', outfmt: '6', args
  } });
  let reply = { keys: ['a'.repeat(64)] };
  const generate = (manifest, args = []) => block({
    useProteinBlastp: true, preparedJobs: [pair(args)], throwIfGenerationCanceled() {},
    useDerivedProteinPayloadCache: false, converterRecordIndexes: [], getSeqEntry() {}, getViewTransform() {},
    async runDiagramHelperOperation(operation, payload) {
      requests.push({ operation, payload });
      return { result: reply };
    },
    DIAGRAM_HELPER_OPERATIONS: { BUILD_PROTEIN_LOSAT_CACHE_KEYS: 'buildProteinLosatCacheKeys' },
    cloneJsonData: structuredClone, workingProteinIdentityManifest: manifest, normalizeLosatArgs: (value) => value,
    derivedSettings: {}, pairNumbers: [], requireProteinLosatCacheKeys, ...probe
  });
  const manifest = mergeProteinIdentityManifests([]);
  assert.deepEqual(await generate(manifest), { proteinCacheKeys: ['a'.repeat(64)], derivedCacheKey: '' });
  // A threshold change, past the identical-Generate skip: the same manifest and pairs.
  assert.deepEqual(await generate(manifest), { proteinCacheKeys: ['a'.repeat(64)], derivedCacheKey: '' });
  assert.equal(requests.length, 1, 'the second Generate reuses the reply');
  assert.equal(requests[0].payload.identityManifest, manifest, 'the frozen manifest is sent without a copy');
  // Every input of the reply is in the key: other options, and another manifest, ask again.
  reply = { keys: ['b'.repeat(64)] };
  assert.deepEqual((await generate(manifest, ['--max-target-seqs', '6'])).proteinCacheKeys, ['b'.repeat(64)]);
  assert.deepEqual((await generate(mergeProteinIdentityManifests([emptyProteinIdentityManifest()]))).proteinCacheKeys,
    ['b'.repeat(64)]);
  assert.equal(requests.length, 3);
  // A refused reply is not kept.
  const other = mergeProteinIdentityManifests([emptyProteinIdentityManifest(), emptyProteinIdentityManifest()]);
  reply = { error: new Error('refused') };
  await assert.rejects(generate(other), /refused/);
  reply = { keys: ['c'.repeat(64)] };
  assert.deepEqual((await generate(other)).proteinCacheKeys, ['c'.repeat(64)]);
  assert.equal(requests.length, 5);
});
