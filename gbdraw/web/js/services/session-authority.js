// @ts-check
import { validateAnnotationWarnings } from './session-feature-metadata.js';
import { validateComparisonWarnings } from './comparison-warnings.js';
import { assertSafeObjectKeys } from './safe-object-keys.js';
import { validateWebFileBindings } from './session-resource-backing.js';
import { validateCurrentWriterActiveConfig, validateAlignmentResetReceiptShape } from './session-active-config-contract.js';
import { migrateLegacyLinearLabelVisibility } from './linear-label-visibility.js';
import { migrateLegacyRecordDisplayDrafts } from './record-display-model.js';
import { canonicalFeatureOverrides, validateFeatureIdentityNotices } from './feature-placement.js';
import { RENDERED_ID_FEATURE_EDIT_FIELDS, migrateSessionFeaturePlacements } from './feature-edit-migration.js';
import { collectCanonicalResourceIds } from './canonical-resource-references.js';

// Session 44 introduced the current active-config and record-display draft
// shapes; Session 45 keys per-feature edits by source identity (design Q4).
export const TYPED_DRAFT_SESSION_VERSION = 44;
export const FEATURE_IDENTITY_SESSION_VERSION = 45;
const FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION = Object.freeze({ 44: 4, 45: 5 });

export const SESSION_TOP_LEVEL_AUTHORITY = Object.freeze({
  format: 'document',
  version: 'document',
  createdAt: 'document',
  title: 'document',
  renderRequest: 'canonical-render',
  resources: 'resource',
  webFiles: 'resource-binding',
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
  otherModeResult: 'artifact',
  cliInvocation: 'provenance'
});

const WEB_EDITOR_UI_FIELDS = Object.freeze([
  'cInputType',
  'lInputType',
  'zoom',
  'canvasPan',
  'canvasPadding',
  'selectedResultIndex',
  'featurePanelTab',
  'downloadDpi',
  'autoLabelReflow',
  'linearTypographyLinked',
  'paletteInstantPreviewEnabled',
  'appliedPaletteName',
  'appliedPaletteColors',
  'pendingPaletteName',
  'pendingPaletteColors',
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

// Older Sessions' rendered-ID edit maps are read once and migrated to
// `featureOverrides` (services/feature-edit-migration.js).
const ARTIFACT_FEATURE_FIELDS = Object.freeze([
  'extractedFeatures',
  'biologicalFeatures',
  'featureRecordIds',
  'selectedFeatureRecordIdx',
  'featureColorOverrides',
  'featureVisibilityManualRules',
  'featureOverrides',
  'labelOverrideRows',
  'labelTextBulkOverrides',
  ...RENDERED_ID_FEATURE_EDIT_FIELDS
]);

const isPlainObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const adoptedCurrentDocuments = new WeakSet();
const adoptedCanonicalOwners = new WeakSet();

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

/**
 * @param {Record<string, any>} sessionData Unvalidated Session data.
 * @returns {void}
 */
export const validateCurrentComparisonAuthority = (sessionData) => {
  const config = isPlainObject(sessionData.config) ? sessionData.config : {};
  const ui = isPlainObject(sessionData.ui) ? sessionData.ui : {};
  const webFiles = isPlainObject(sessionData.webFiles) ? sessionData.webFiles : {};
  const bindings = webFiles.bindings || {};

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

const hasInput = value => Array.isArray(value) ? value.some(hasInput) : value != null;

// The same complete inventory is used before Save and at document admission.
// An inactive input remains biological even when it cannot render the active mode.
/**
 * @param {Record<string, any>} [files] The Session's Web input inventory.
 * @returns {boolean}
 */
export const hasBiologicalSessionInputs = (files = {}) => (
  ['c_gb', 'c_gff', 'c_fasta', 'c_conservation_fastas', 'c_conservation_sequence_sources']
    .some(key => hasInput(files[key]))
  || (files.linearSeqs || []).some(row => ['gb', 'gff', 'fasta'].some(key => hasInput(row[key])))
);

/**
 * @param {Record<string, any> | null | undefined} data
 * @returns {boolean}
 */
export const isSettingsOnlySessionDocument = data => data != null && [42, 44, 45].includes(data.version)
  && Object.hasOwn(data, 'renderRequest') && data.renderRequest === null;

const validateSettingsOnlyDocument = data => {
  const bindings = data.webFiles?.bindings;
  if (!isPlainObject(data.resources) || !isPlainObject(bindings)
    || !Array.isArray(bindings.linearSeqs)
    || ['c_gb', 'c_gff', 'c_fasta'].some(key => !Object.hasOwn(bindings, key))) {
    throw new Error('Settings-only Session requires resources and an explicit Web input inventory.');
  }
  if (hasBiologicalSessionInputs(bindings)
    || Object.keys(data.webFiles).some(key => key !== 'bindings')) {
    throw new Error('Settings-only Session cannot contain biological sources.');
  }
  if (data.results?.length !== 0 || data.editorState?.featureCatalog !== null
    || data.cliInvocation != null || Object.keys(data.runMetadata || {}).length
    || (data.losatCache?.entries || []).length || (data.losatDerivedCache?.entries || []).length
    || Object.keys(data.legacyArtifacts || {}).length
    || ['proteinSets', 'recordAnalyses', 'recordInstances']
      .some(key => Object.keys(data.proteinIdentityManifest?.[key] || {}).length)
    || bindings.c_conservation_blasts_source === 'losat-cache') {
    throw new Error('Settings-only Session cannot contain committed render artifacts.');
  }
  let storedConfig = data.version >= TYPED_DRAFT_SESSION_VERSION || !isPlainObject(data.config?.adv)
    ? data.config
    : {
        ...data.config,
        adv: migrateLegacyLinearLabelVisibility(data.config.adv)
      };
  if (data.version < TYPED_DRAFT_SESSION_VERSION && isPlainObject(storedConfig)
    && Object.prototype.hasOwnProperty.call(storedConfig, 'recordDisplayDrafts')) {
    storedConfig = {
      ...storedConfig,
      recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(
        storedConfig.recordDisplayDrafts
      )
    };
  }
  if (data.version < FEATURE_IDENTITY_SESSION_VERSION && isPlainObject(storedConfig)
    && Object.prototype.hasOwnProperty.call(storedConfig, 'featurePlacementOverrides')) {
    storedConfig = {
      ...storedConfig,
      featurePlacementOverrides: migrateSessionFeaturePlacements(storedConfig.featurePlacementOverrides)
    };
  }
  validateCurrentWriterActiveConfig({ mode: data.ui?.mode, storedConfig });
  const referenced = new Set();
  const visit = value => {
    if (!value || typeof value !== 'object') return;
    if (Object.hasOwn(value, 'resourceId')) referenced.add(value.resourceId);
    Object.values(value).forEach(visit);
  };
  visit(bindings);
  if (Object.keys(data.resources).some(id => !referenced.has(id))) {
    throw new Error('Settings-only Session contains an unbound resource.');
  }
};

// One committed Result set: the top-level fields, or `otherModeResult`, whose
// fields mirror them.
/**
 * @param {Record<string, any>} set
 * @param {number | string} version
 */
const validateCommittedResultSet = (set, version) => {
  if (!Array.isArray(set.results)) {
    throw new Error(`Session version ${String(version)} requires a results array.`);
  }
  set.results.forEach((result) => {
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
  const editorState = set.editorState;
  if (
    !isPlainObject(editorState)
    || !Object.prototype.hasOwnProperty.call(editorState, 'featureCatalog')
  ) {
    throw new Error(
      `Session version ${String(version)} requires editorState.featureCatalog.`
    );
  }
  if (set.renderRequest?.schema >= 8
    && set.renderRequest.layout?.similarityAlignment
    && !Object.hasOwn(editorState, 'alignmentResetReceipt')) {
    throw new Error('Current alignment Session requires editorState.alignmentResetReceipt.');
  }
  validateAlignmentResetReceiptShape(editorState.alignmentResetReceipt, set.renderRequest);
  const featureCatalog = editorState.featureCatalog;
  const expectedCatalogSchema = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION[Number(version)] ?? 3;
  if (
    featureCatalog !== null
    && (!isPlainObject(featureCatalog) || featureCatalog.schema !== expectedCatalogSchema)
  ) {
    throw new Error(
      `Session version ${String(version)} requires a schema-${expectedCatalogSchema} editorState.featureCatalog.`
    );
  }
  if (
    Array.isArray(set.results)
    && set.results.length > 0
    && featureCatalog === null
  ) {
    throw new Error(
      `Session version ${String(version)} requires a feature catalog for saved results.`
    );
  }
};

const OTHER_MODE_RESULT_FIELDS = new Set([
  'renderRequest', 'results', 'editorState', 'ui', 'runMetadata', 'cliInvocation'
]);
// The per-set part of the shared `ui` and `editorState` objects.
const OTHER_MODE_RESULT_UI_FIELDS = new Set([
  'selectedResultIndex', 'generatedLegendPosition', 'generatedMultiRecordCanvas',
  'generatedCircularPlotTitlePosition', 'appliedPaletteName', 'appliedPaletteColors'
]);
const OTHER_MODE_RESULT_EDITOR_FIELDS = new Set([
  'featureCatalog', 'alignmentResetReceipt', 'legend', 'originalSvgStroke'
]);
const OTHER_MODE_RESULT_LEGEND_FIELDS = new Set(['originalOrder', 'originalColors']);
const OTHER_MODE_RESULT_RUN_METADATA_FIELDS = new Set([
  'trackSlotGeometry', 'annotationWarnings', 'featureIdentityNotices', 'comparisonWarnings'
]);

// E1: Session 45 keeps the other diagram mode's Result set beside the
// top-level set. It needs a committed top-level request of the other mode, at
// least one Result, and resources in the top-level table.
/**
 * @param {Record<string, any>} sessionData
 * @param {number | string} version
 */
const validateOtherModeResult = (sessionData, version) => {
  const other = sessionData.otherModeResult;
  if (!isPlainObject(other) || Object.keys(other).some((key) => !OTHER_MODE_RESULT_FIELDS.has(key))) {
    throw new Error('Session otherModeResult must contain only a committed Result set.');
  }
  const mode = other.renderRequest?.mode;
  if (!isPlainObject(sessionData.renderRequest) || !isPlainObject(other.renderRequest)
    || !['circular', 'linear'].includes(mode) || mode === sessionData.renderRequest.mode
    || other.renderRequest.schema !== sessionData.renderRequest.schema) {
    throw new Error('Session otherModeResult requires a committed request of the other mode.');
  }
  validateCommittedResultSet(other, version);
  if (other.results.length === 0) throw new Error('Session otherModeResult requires a Result.');
  const missing = [...collectCanonicalResourceIds(other.renderRequest)]
    .filter((resourceId) => !Object.hasOwn(sessionData.resources || {}, resourceId));
  if (missing.length > 0) {
    throw new Error(`Session otherModeResult names missing resource(s): ${missing.join(', ')}`);
  }
  // One rule with gbdraw/session_io.py: a present field is an object (null is not).
  const ui = Object.hasOwn(other, 'ui') ? other.ui : {};
  const legend = Object.hasOwn(other.editorState, 'legend') ? other.editorState.legend : {};
  const runMetadata = Object.hasOwn(other, 'runMetadata') ? other.runMetadata : {};
  if (!isPlainObject(ui) || Object.keys(ui).some((key) => !OTHER_MODE_RESULT_UI_FIELDS.has(key))
    || Object.keys(other.editorState).some((key) => !OTHER_MODE_RESULT_EDITOR_FIELDS.has(key))
    || !isPlainObject(legend) || Object.keys(legend).some((key) => !OTHER_MODE_RESULT_LEGEND_FIELDS.has(key))) {
    throw new Error('Session otherModeResult editorState and ui hold only that Result set\'s fields.');
  }
  if (!isPlainObject(runMetadata)
    || Object.keys(runMetadata).some((key) => !OTHER_MODE_RESULT_RUN_METADATA_FIELDS.has(key))) {
    throw new Error('Session otherModeResult.runMetadata holds only that Result set\'s metadata.');
  }
  validateAnnotationWarnings(runMetadata.annotationWarnings, other.results);
  validateFeatureIdentityNotices(runMetadata.featureIdentityNotices, other.results);
  validateComparisonWarnings(runMetadata.comparisonWarnings, other.results);
};

/**
 * @param {Record<string, any>} sessionData Unvalidated Session data.
 * @param {number | string} version
 * @returns {void}
 */
export const validateSessionAuthorityInventory = (sessionData, version) => {
  if (!sessionData || typeof sessionData !== 'object' || Array.isArray(sessionData)) {
    throw new Error('Session authority inventory requires an object.');
  }
  assertSafeObjectKeys(sessionData, 'Session');
  if (Object.hasOwn(sessionData, 'otherModeResult') && Number(version) < FEATURE_IDENTITY_SESSION_VERSION) {
    throw new Error(`Session version ${String(version)} cannot contain otherModeResult.`);
  }
  const bindings = validateWebFileBindings(sessionData.webFiles, sessionData.resources);
  if (bindings?.schema === 2 && ![41, 42, 44, 45].includes(Number(version))) {
    throw new Error('Web binding schema 2 requires session version 41, 42, 44, or 45.');
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
    for (const field of ['cInputType', 'lInputType']) {
      if (
        Object.prototype.hasOwnProperty.call(sessionData.ui || {}, field)
        && !['gb', 'gff'].includes(sessionData.ui[field])
      ) {
        throw new Error(`Session ui.${field} must be gb or gff when present.`);
      }
    }
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
    if (Number(version) >= FEATURE_IDENTITY_SESSION_VERSION) {
      const retired = RENDERED_ID_FEATURE_EDIT_FIELDS.filter((field) => Object.hasOwn(features, field));
      if (retired.length > 0) {
        throw new Error(
          `Session version ${String(version)} cannot contain rendered-ID feature edits: ${retired.join(', ')}`
        );
      }
      if (Object.hasOwn(features, 'featureOverrides')) {
        if (!isPlainObject(features.featureOverrides)) {
          throw new Error('Session features.featureOverrides must be an object.');
        }
        canonicalFeatureOverrides(features.featureOverrides);
      }
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
    validateCommittedResultSet(sessionData, version);
  }
  if (Object.hasOwn(sessionData, 'otherModeResult')) validateOtherModeResult(sessionData, version);
  const unknown = Object.keys(sessionData).filter(
    (key) => !Object.prototype.hasOwnProperty.call(SESSION_TOP_LEVEL_AUTHORITY, key)
  );
  if (unknown.length > 0) {
    throw new Error(`Session contains unclassified top-level field(s): ${unknown.join(', ')}`);
  }
  validateAnnotationWarnings(sessionData.runMetadata?.annotationWarnings, sessionData.results);
  validateFeatureIdentityNotices(sessionData.runMetadata?.featureIdentityNotices, sessionData.results);
  validateComparisonWarnings(sessionData.runMetadata?.comparisonWarnings, sessionData.results);
  if (isSettingsOnlySessionDocument(sessionData)) validateSettingsOnlyDocument(sessionData);
};

/**
 * @param {Record<string, any>} canonical `{ renderRequest, resources, webFiles }`.
 * @returns {Record<string, any>} The same object, adopted.
 */
export const adoptRuntimeCanonicalSession = (canonical) => {
  if (
    !isPlainObject(canonical)
    || !isPlainObject(canonical.renderRequest)
    || !isPlainObject(canonical.resources)
    || (
      canonical.webFiles !== undefined
      && !isPlainObject(canonical.webFiles)
    )
  ) {
    throw new Error('A canonical render request is required for adoptive ownership.');
  }
  adoptedCanonicalOwners.add(canonical);
  return canonical;
};

/**
 * @param {Record<string, any>} sessionData Unvalidated Session data.
 * @param {number} currentVersion
 * @returns {{ document: Record<string, any>, canonical: Record<string, any> | null }}
 */
export const adoptCurrentSessionDocument = (sessionData, currentVersion) => {
  validateSessionAuthorityInventory(sessionData, currentVersion);
  if (sessionData.version !== currentVersion) {
    throw new Error('Only the current session schema can use adoptive ownership.');
  }
  const canonical = isSettingsOnlySessionDocument(sessionData) ? null : adoptRuntimeCanonicalSession({
    renderRequest: sessionData.renderRequest,
    resources: sessionData.resources,
    webFiles: isPlainObject(sessionData.webFiles) ? sessionData.webFiles : {}
  });
  adoptedCurrentDocuments.add(sessionData);
  return { document: sessionData, canonical };
};

/**
 * @param {any} value
 * @returns {boolean}
 */
export const isAdoptedCurrentSessionDocument = (value) => (
  Boolean(value) && adoptedCurrentDocuments.has(value)
);

/**
 * @param {any} value
 * @returns {boolean}
 */
export const isAdoptedCanonicalSession = (value) => (
  Boolean(value) && adoptedCanonicalOwners.has(value)
);

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {{ ui: Record<string, any> }}
 */
export const projectWebOnlyEditorMetadata = (sessionData) => ({
  ui: copyFields(sessionData?.ui, WEB_EDITOR_UI_FIELDS)
});

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {Record<string, any>}
 */
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

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {{ format: any, version: any, createdAt: any, title: string }}
 */
export const projectDocumentMetadata = (sessionData) => ({
  format: sessionData?.format,
  version: sessionData?.version,
  createdAt: sessionData?.createdAt,
  title: typeof sessionData?.title === 'string' ? sessionData.title : ''
});
