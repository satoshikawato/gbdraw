import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-session-authority-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempRoot, 'session-authority.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'session-authority.js'), 'utf8'),
  'utf8'
);

const {
  SESSION_TOP_LEVEL_AUTHORITY,
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
  validateSessionAuthorityInventory
} = await import(pathToFileURL(join(tempRoot, 'session-authority.js')));

assert.deepEqual(Object.keys(SESSION_TOP_LEVEL_AUTHORITY).sort(), [
  'cliInvocation', 'config', 'createdAt', 'editorState', 'features', 'files', 'format',
  'legacyArtifacts', 'losatCache', 'losatDerivedCache', 'orthogroupState',
  'proteinIdentityManifest', 'renderRequest', 'resources',
  'results', 'title', 'ui', 'version', 'webFiles'
].sort());

const session = {
  format: 'gbdraw-session', version: 33, createdAt: 'now', title: 'Canonical',
  renderRequest: {}, resources: {}, webFiles: {}, config: {}, files: {},
  ui: {
    mode: 'linear', legend: 'left', linearPlotTitlePosition: 'top', zoom: 1.5,
    canvasPan: { x: 3, y: 4 }, generatedLegendPosition: 'right', downloadDpi: 300
  },
  features: {
    extractedFeatures: [{ id: 'f1' }], featureColorOverrides: { f1: '#ffffff' },
    featureVisibilityManualRules: [{ action: 'off' }], labelTextFeatureOverrides: { f1: 'stored' }
  },
  editorState: { legend: { entries: [] } }, results: [{ name: 'preview', content: '<svg/>' }],
  orthogroupState: {}, losatCache: {}, losatDerivedCache: {},
  proteinIdentityManifest: { schema: 1 },
  legacyArtifacts: { proteinRawCandidates: { schema: 1, entries: [] } },
  cliInvocation: null
};
validateSessionAuthorityInventory(session, 33);
assert.deepEqual(projectWebOnlyEditorMetadata(session).ui, {
  zoom: 1.5,
  canvasPan: { x: 3, y: 4 },
  downloadDpi: 300
});
assert.deepEqual(projectArtifactState(session).features, {
  extractedFeatures: [{ id: 'f1' }]
});
assert.deepEqual(projectArtifactState(session).ui, { generatedLegendPosition: 'right' });
assert.deepEqual(projectArtifactState(session).proteinIdentityManifest, { schema: 1 });
assert.deepEqual(projectArtifactState(session).legacyArtifacts, session.legacyArtifacts);
assert.equal(projectDocumentMetadata(session).title, 'Canonical');
assert.throws(
  () => validateSessionAuthorityInventory({ ...session, unknownField: true }, 33),
  /unclassified top-level field.*unknownField/
);
assert.doesNotThrow(() => validateSessionAuthorityInventory({ ...session, unknownField: true }, 30));

