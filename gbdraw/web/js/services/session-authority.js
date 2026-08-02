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

const LINEAR_COMPARISON_PLAN_MODES = new Set(['none', 'adjacent', 'selected']);
const LINEAR_COMPARISON_SOURCES = new Set(['losat', 'upload']);
const LINEAR_COMPARISON_EDGE_FIELDS = new Set([
  'id',
  'queryUid',
  'subjectUid',
  'included',
  'fileActive',
  'losatFilenameActive',
  'source',
  'losatFilename'
]);

const assertNoOwnField = (value, field, message) => {
  if (isPlainObject(value) && Object.prototype.hasOwnProperty.call(value, field)) {
    throw new Error(message);
  }
};

const validateLinearComparisonPlan = (plan) => {
  if (!isPlainObject(plan)) {
    throw new Error('Current Web comparison draft requires config.linearComparisonPlan.');
  }
  if (!LINEAR_COMPARISON_PLAN_MODES.has(plan.mode)) {
    throw new Error('config.linearComparisonPlan.mode is invalid.');
  }
  if (!LINEAR_COMPARISON_SOURCES.has(plan.defaultSource)) {
    throw new Error('config.linearComparisonPlan.defaultSource is invalid.');
  }
  if (!Array.isArray(plan.edges)) {
    throw new Error('config.linearComparisonPlan.edges must be an array.');
  }
  const ids = new Set();
  plan.edges.forEach((edge) => {
    if (!isPlainObject(edge)) {
      throw new Error('Each config.linearComparisonPlan edge must be an object.');
    }
    const unknown = Object.keys(edge).filter((field) => !LINEAR_COMPARISON_EDGE_FIELDS.has(field));
    if (unknown.length > 0) {
      throw new Error(
        `config.linearComparisonPlan edge contains retired or unknown field(s): ${unknown.join(', ')}`
      );
    }
    const id = typeof edge.id === 'string' ? edge.id.trim() : '';
    const queryUid = typeof edge.queryUid === 'string' ? edge.queryUid.trim() : '';
    const subjectUid = typeof edge.subjectUid === 'string' ? edge.subjectUid.trim() : '';
    if (!id || !queryUid || !subjectUid) {
      throw new Error('Current comparison-plan edges require stable IDs and endpoint UIDs.');
    }
    if (ids.has(id)) {
      throw new Error(`Current comparison-plan edge ID is duplicated: ${id}.`);
    }
    ids.add(id);
    if (
      typeof edge.included !== 'boolean'
      || typeof edge.fileActive !== 'boolean'
      || typeof edge.losatFilenameActive !== 'boolean'
      || !LINEAR_COMPARISON_SOURCES.has(edge.source)
      || typeof edge.losatFilename !== 'string'
    ) {
      throw new Error('Current comparison-plan edge metadata is invalid.');
    }
  });
  return ids;
};

export const validateCurrentComparisonAuthority = (sessionData) => {
  const config = isPlainObject(sessionData.config) ? sessionData.config : {};
  const ui = isPlainObject(sessionData.ui) ? sessionData.ui : {};
  const webFiles = isPlainObject(sessionData.webFiles) ? sessionData.webFiles : {};
  const hasBindings = Object.prototype.hasOwnProperty.call(webFiles, 'bindings');
  if (hasBindings && !isPlainObject(webFiles.bindings)) {
    throw new Error('Session webFiles.bindings must be an object.');
  }
  const bindings = hasBindings ? webFiles.bindings : {};
  if (hasBindings && bindings.schema !== 1) {
    throw new Error('Unsupported Web file binding schema.');
  }

  assertNoOwnField(config, 'blastSource', 'Current sessions cannot contain config.blastSource.');
  assertNoOwnField(config.adv, 'blastSource', 'Current sessions cannot contain config.adv.blastSource.');
  assertNoOwnField(ui, 'blastSource', 'Current sessions cannot contain ui.blastSource.');
  assertNoOwnField(
    config.linearRecordLayout,
    'comparisons',
    'Current sessions cannot contain config.linearRecordLayout.comparisons.'
  );

  (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach((sequence) => {
    assertNoOwnField(
      sequence,
      'blast',
      'Current sessions cannot contain per-record BLAST file bindings.'
    );
    assertNoOwnField(
      sequence,
      'losat_filename',
      'Current sessions cannot contain per-record LOSAT filenames.'
    );
  });
  (Array.isArray(webFiles.linearRecordMetadata) ? webFiles.linearRecordMetadata : [])
    .forEach((metadata) => {
      assertNoOwnField(
        metadata,
        'losatFilename',
        'Current sessions cannot contain per-record LOSAT filename metadata.'
      );
      assertNoOwnField(
        metadata,
        'losat_filename',
        'Current sessions cannot contain per-record LOSAT filename metadata.'
      );
    });
  assertNoOwnField(
    bindings,
    'linearCanonicalComparisons',
    'Current sessions cannot bind comparison artifacts outside the committed request.'
  );

  const hasWebComparisonDraft = Object.prototype.hasOwnProperty.call(
    config,
    'linearRecordLayout'
  ) || Object.prototype.hasOwnProperty.call(config, 'linearComparisonPlan');
  if (
    Object.prototype.hasOwnProperty.call(bindings, 'linearComparisons')
    && !Array.isArray(bindings.linearComparisons)
  ) {
    throw new Error('Current comparison file bindings must be an array.');
  }
  const bindingEntries = Array.isArray(bindings.linearComparisons) ? bindings.linearComparisons : [];
  if (!hasWebComparisonDraft && bindingEntries.length === 0) return;

  const edgeIds = validateLinearComparisonPlan(config.linearComparisonPlan);
  const resources = isPlainObject(sessionData.resources) ? sessionData.resources : {};
  const boundIds = new Set();
  bindingEntries.forEach((binding) => {
    if (!isPlainObject(binding)) {
      throw new Error('Each current comparison file binding must be an object.');
    }
    const unknown = Object.keys(binding).filter((field) => !['id', 'file'].includes(field));
    if (unknown.length > 0) {
      throw new Error(
        `Current comparison file binding duplicates plan metadata: ${unknown.join(', ')}`
      );
    }
    const id = typeof binding.id === 'string' ? binding.id.trim() : '';
    if (!id || !edgeIds.has(id)) {
      throw new Error('Current comparison file binding must reference a plan edge ID.');
    }
    if (boundIds.has(id)) {
      throw new Error(`Current comparison file binding is duplicated: ${id}.`);
    }
    if (!isPlainObject(binding.file)) {
      throw new Error('Each current comparison file binding requires a file resource binding.');
    }
    const resourceId = typeof binding.file.resourceId === 'string'
      ? binding.file.resourceId.trim()
      : '';
    if (!resourceId) {
      throw new Error('Each current comparison file binding requires a resourceId.');
    }
    if (!Object.prototype.hasOwnProperty.call(resources, resourceId)) {
      throw new Error(`Current comparison file binding references a missing resource: ${resourceId}.`);
    }
    boundIds.add(id);
  });
  config.linearComparisonPlan.edges.forEach((edge) => {
    if (edge.fileActive && !boundIds.has(edge.id)) {
      throw new Error(
        `Active comparison file is missing its Web file binding: ${edge.id}.`
      );
    }
  });
};

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
    validateCurrentComparisonAuthority(sessionData);
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
