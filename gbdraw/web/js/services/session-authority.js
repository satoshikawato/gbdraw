export const SESSION_TOP_LEVEL_AUTHORITY = Object.freeze({
  format: 'document',
  version: 'document',
  createdAt: 'document',
  title: 'document',
  renderRequest: 'canonical-render',
  resources: 'canonical-render',
  webFiles: 'canonical-render',
  config: 'legacy-or-editor-metadata',
  ui: 'editor-metadata',
  files: 'legacy-fallback',
  results: 'artifact',
  features: 'artifact-or-legacy-semantic',
  editorState: 'artifact',
  orthogroupState: 'artifact',
  losatCache: 'artifact',
  losatDerivedCache: 'artifact',
  proteinIdentityManifest: 'artifact',
  legacyArtifacts: 'artifact',
  cliInvocation: 'provenance'
});

const WEB_EDITOR_UI_FIELDS = Object.freeze([
  'zoom',
  'canvasPan',
  'canvasPadding',
  'selectedResultIndex',
  'featurePanelTab',
  'downloadDpi',
  'autoLabelReflow',
  'paletteInstantPreviewEnabled'
]);

const ARTIFACT_UI_FIELDS = Object.freeze([
  'generatedLegendPosition',
  'generatedMultiRecordCanvas',
  'generatedCircularPlotTitlePosition'
]);

const ARTIFACT_FEATURE_FIELDS = Object.freeze([
  'extractedFeatures',
  'featureSelectorSafetyScope',
  'featureRecordIds',
  'selectedFeatureRecordIdx'
]);

const copyFields = (source, fields) => {
  const projected = {};
  if (!source || typeof source !== 'object' || Array.isArray(source)) return projected;
  fields.forEach((field) => {
    if (Object.prototype.hasOwnProperty.call(source, field)) {
      projected[field] = source[field];
    }
  });
  return projected;
};

export const validateSessionAuthorityInventory = (sessionData, version) => {
  if (!sessionData || typeof sessionData !== 'object' || Array.isArray(sessionData)) {
    throw new Error('Session authority inventory requires an object.');
  }
  if (Number(version) < 31) return;
  const unknown = Object.keys(sessionData).filter(
    (key) => !Object.prototype.hasOwnProperty.call(SESSION_TOP_LEVEL_AUTHORITY, key)
  );
  if (unknown.length > 0) {
    throw new Error(`Session contains unclassified top-level field(s): ${unknown.join(', ')}`);
  }
};

export const projectWebOnlyEditorMetadata = (sessionData) => ({
  ui: copyFields(sessionData?.ui, WEB_EDITOR_UI_FIELDS)
});

export const projectArtifactState = (sessionData) => ({
  results: Array.isArray(sessionData?.results) ? sessionData.results : [],
  ui: copyFields(sessionData?.ui, ARTIFACT_UI_FIELDS),
  features: copyFields(sessionData?.features, ARTIFACT_FEATURE_FIELDS),
  editorState: sessionData?.editorState || {},
  orthogroupState: sessionData?.orthogroupState || {},
  losatCache: sessionData?.losatCache || {},
  losatDerivedCache: sessionData?.losatDerivedCache || {},
  proteinIdentityManifest: sessionData?.proteinIdentityManifest || {},
  legacyArtifacts: sessionData?.legacyArtifacts || {}
});

export const projectDocumentMetadata = (sessionData) => ({
  format: sessionData?.format,
  version: sessionData?.version,
  createdAt: sessionData?.createdAt,
  title: typeof sessionData?.title === 'string' ? sessionData.title : ''
});

