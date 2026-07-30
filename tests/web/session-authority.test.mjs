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
  'results', 'runMetadata', 'title', 'ui', 'version', 'webFiles'
].sort());

const session = {
  format: 'gbdraw-session', version: 39, createdAt: 'now', title: 'Canonical',
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
  proteinIdentityManifest: {
    schema: 2,
    proteinSets: {},
    recordAnalyses: {},
    recordInstances: {}
  },
  legacyArtifacts: { proteinRawCandidates: { schema: 1, entries: [] } },
  runMetadata: { trackSlotGeometry: { schema: 1, records: [] } },
  cliInvocation: null
};
validateSessionAuthorityInventory(session, 39);
const currentSession = {
  ...session,
  version: 40,
  features: {
    featureColorOverrides: { f1: '#ffffff' },
    featureVisibilityManualRules: [{ action: 'off' }],
    labelTextFeatureOverrides: { f1: 'stored' }
  },
  editorState: {
    legend: { entries: [] },
    featureCatalog: { schema: 3, items: [] }
  }
};
delete currentSession.files;
assert.doesNotThrow(() => validateSessionAuthorityInventory(currentSession, 40));
assert.doesNotThrow(() => validateSessionAuthorityInventory({
  ...currentSession,
  results: [],
  editorState: { featureCatalog: null }
}, 40));
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: 'not-an-array',
    editorState: { featureCatalog: null }
  }, 40),
  /requires a results array/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: [],
    editorState: {}
  }, 40),
  /requires editorState.featureCatalog/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: []
  }, 40),
  /Session features must be an object/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    orthogroupState: []
  }, 40),
  /Session orthogroupState must be an object/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: [{
      name: 'preview.interactive.svg',
      content: '<svg><metadata id="gbdraw-interactive-feature-metadata"/></svg>'
    }],
    editorState: {
      featureCatalog: { schema: 3, items: [] }
    }
  }, 40),
  /named plain SVG/
);
assert.throws(
  () => validateSessionAuthorityInventory({ ...currentSession, files: {} }, 40),
  /version 40 cannot contain legacy files/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: {
      ...currentSession.features,
      featureCatalog: {
        schema: 1,
        encoding: 'biological-authority-v1',
        extracted: []
      }
    }
  }, 40),
  /branch-only feature field.*featureCatalog/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: {
      ...currentSession.features,
      biologicalFeatures: []
    }
  }, 40),
  /branch-only feature field.*biologicalFeatures/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    orthogroupState: { groups: [] }
  }, 40),
  /cannot contain duplicated orthogroup groups/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    editorState: { featureCatalog: { schema: 1 } }
  }, 40),
  /requires a schema-3 editorState.featureCatalog/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    editorState: { featureCatalog: null }
  }, 40),
  /requires a feature catalog for saved results/
);
assert.deepEqual(projectWebOnlyEditorMetadata(session).ui, {
  legend: 'left',
  linearPlotTitlePosition: 'top',
  zoom: 1.5,
  canvasPan: { x: 3, y: 4 },
  downloadDpi: 300
});
const layoutPreferences = {
  circular: {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: 'right', plotTitlePosition: 'bottom' }
  },
  linear: { legend: 'bottom', plotTitlePosition: 'top' }
};
assert.deepEqual(
  projectWebOnlyEditorMetadata({ ui: { layoutPreferences } }).ui,
  { layoutPreferences }
);
assert.deepEqual(projectArtifactState(session).features, {
  extractedFeatures: [{ id: 'f1' }],
  featureColorOverrides: { f1: '#ffffff' },
  featureVisibilityManualRules: [{ action: 'off' }],
  labelTextFeatureOverrides: { f1: 'stored' }
});
assert.deepEqual(projectArtifactState(session).ui, { generatedLegendPosition: 'right' });
assert.deepEqual(
  projectArtifactState(session).proteinIdentityManifest,
  session.proteinIdentityManifest
);
assert.deepEqual(projectArtifactState(session).legacyArtifacts, session.legacyArtifacts);
assert.deepEqual(projectArtifactState(session).runMetadata, session.runMetadata);
assert.equal(projectDocumentMetadata(session).title, 'Canonical');
assert.throws(
  () => validateSessionAuthorityInventory({ ...session, unknownField: true }, 39),
  /unclassified top-level field.*unknownField/
);
assert.throws(
  () => validateSessionAuthorityInventory({ ...currentSession, unknownField: true }, 40),
  /unclassified top-level field.*unknownField/
);
assert.doesNotThrow(() => validateSessionAuthorityInventory({ ...session, unknownField: true }, 30));
