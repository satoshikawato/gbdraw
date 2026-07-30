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
  runMetadata: 'artifact',
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
  'paletteInstantPreviewEnabled',
  'layoutPreferences',
  // Older editor payloads stored the same preference state in parallel fields.
  'legend',
  'circularLegendPosition',
  'linearLegendPosition',
  'circularPlotTitlePosition',
  'linearPlotTitlePosition',
  'circularSingleRecordLegendPosition',
  'circularSingleRecordPlotTitlePosition',
  'circularMultiRecordLegendPosition',
  'circularMultiRecordPlotTitlePosition'
]);

const ARTIFACT_UI_FIELDS = Object.freeze([
  'generatedLegendPosition',
  'generatedMultiRecordCanvas',
  'generatedCircularPlotTitlePosition'
]);

const ARTIFACT_FEATURE_FIELDS = Object.freeze([
  'extractedFeatures',
  'biologicalFeatures',
  'featureSelectorSafetyScope',
  'featureRecordIds',
  'selectedFeatureRecordIdx',
  'featureColorOverrides',
  'featureVisibilityManualRules',
  'featureVisibilityOverrides',
  'labelTextFeatureOverrides',
  'labelOverrideRows',
  'labelTextBulkOverrides',
  'labelTextFeatureOverrideSources',
  'labelVisibilityOverrides',
  'labelOverrideContextKey'
]);

const isPlainObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS = Object.freeze([
  'extractedFeatures',
  'biologicalFeatures',
  'featureSelectorSafetyScope',
  'featureRecordIds',
  'featureCatalog'
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
  if (
    Number(version) >= 40 &&
    Object.prototype.hasOwnProperty.call(sessionData, 'files')
  ) {
    throw new Error(
      `Session version ${String(version)} cannot contain legacy files; use resources and webFiles.`
    );
  }
  if (Number(version) >= 40) {
    if (
      Object.prototype.hasOwnProperty.call(sessionData, 'features')
      && !isPlainObject(sessionData.features)
    ) {
      throw new Error('Session features must be an object when present.');
    }
    const features = isPlainObject(sessionData.features) ? sessionData.features : {};
    const forbiddenFeatureFields = CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS.filter(
      (field) => Object.prototype.hasOwnProperty.call(features, field)
    );
    if (forbiddenFeatureFields.length > 0) {
      throw new Error(
        `Session version ${String(version)} contains branch-only feature field(s): `
        + forbiddenFeatureFields.join(', ')
      );
    }
    if (
      Object.prototype.hasOwnProperty.call(sessionData, 'orthogroupState')
      && !isPlainObject(sessionData.orthogroupState)
    ) {
      throw new Error('Session orthogroupState must be an object when present.');
    }
    if (
      isPlainObject(sessionData.orthogroupState)
      && Object.prototype.hasOwnProperty.call(sessionData.orthogroupState, 'groups')
    ) {
      throw new Error(
        `Session version ${String(version)} cannot contain duplicated orthogroup groups.`
      );
    }
    if (!Array.isArray(sessionData.results)) {
      throw new Error(`Session version ${String(version)} requires a results array.`);
    }
    sessionData.results.forEach((result) => {
      const resultName = (
        isPlainObject(result) && typeof result.name === 'string'
          ? result.name.trim()
          : ''
      );
      const content = isPlainObject(result) ? result.content : null;
      if (
        !resultName
        || resultName.toLowerCase().endsWith('.interactive.svg')
        || typeof content !== 'string'
        || !content.includes('<svg')
        || content.includes('gbdraw-interactive-feature-metadata')
        || content.includes('gbdraw-interactive-feature-script')
      ) {
        throw new Error('Each current Session Result must be a named plain SVG.');
      }
    });
    const editorState = sessionData.editorState;
    if (
      !isPlainObject(editorState)
      || !Object.prototype.hasOwnProperty.call(editorState, 'featureCatalog')
    ) {
      throw new Error(
        `Session version ${String(version)} requires editorState.featureCatalog.`
      );
    }
    const featureCatalog = editorState.featureCatalog;
    if (
      featureCatalog !== null
      && (!isPlainObject(featureCatalog) || featureCatalog.schema !== 3)
    ) {
      throw new Error(
        `Session version ${String(version)} requires a schema-3 editorState.featureCatalog.`
      );
    }
    if (
      Array.isArray(sessionData.results)
      && sessionData.results.length > 0
      && featureCatalog === null
    ) {
      throw new Error(
        `Session version ${String(version)} requires a feature catalog for saved results.`
      );
    }
  }
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
  legacyArtifacts: sessionData?.legacyArtifacts || {},
  runMetadata: sessionData?.runMetadata || {}
});

export const projectDocumentMetadata = (sessionData) => ({
  format: sessionData?.format,
  version: sessionData?.version,
  createdAt: sessionData?.createdAt,
  title: typeof sessionData?.title === 'string' ? sessionData.title : ''
});
