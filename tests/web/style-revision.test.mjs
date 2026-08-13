import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-style-revision-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const {
  advanceDocumentEpoch,
  advanceStyleRevision,
  captureStyleSnapshot,
  styleFingerprint,
  styleSnapshotIsCurrent
} = await import(pathToFileURL(join(tempRoot, 'js', 'services', 'style-revision.js')));

const ref = (value) => ({ value });
const state = {
  manualSpecificRules: [{ feat: 'CDS', qual: 'gene', val: '^a$', color: '#ABCDEF', cap: 'A' }],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ tRNA: '#222222', CDS: '#111111' }),
  documentEpoch: ref(4),
  resultGenerationKey: ref(7),
  semanticStyleRevision: ref(2),
  semanticStyleFingerprint: ref(''),
  validatedStyleFingerprintByResultKey: ref(Object.freeze({}))
};

const before = captureStyleSnapshot(state);
assert.equal(before.fingerprint, styleFingerprint({
  rules: state.manualSpecificRules,
  appliedPaletteName: 'default',
  appliedPaletteColors: { CDS: '#111111', tRNA: '#222222' }
}));
assert.equal(styleSnapshotIsCurrent(state, before), true);

state.manualSpecificRules[0].color = '#abcdef';
assert.equal(styleSnapshotIsCurrent(state, before), true, 'hex case is not a semantic change');
state.manualSpecificRules[0].color = '#000000';
assert.equal(styleSnapshotIsCurrent(state, before), false);
assert.throws(
  () => advanceStyleRevision(state, { expected: before }),
  /style changed/
);

state.manualSpecificRules[0].color = '#abcdef';
const transition = advanceStyleRevision(state, {
  expected: before,
  resultFingerprints: { 'result-a': before.fingerprint }
});
assert.equal(transition.revision, 3);
assert.deepEqual(state.validatedStyleFingerprintByResultKey.value, {
  'result-a': before.fingerprint
});

assert.equal(advanceDocumentEpoch(state), 5);
assert.equal(state.semanticStyleRevision.value, 0);
assert.deepEqual(state.validatedStyleFingerprintByResultKey.value, {});
