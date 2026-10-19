// @ts-check
import { diagnosticError, normalizeCaughtError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { describeLoadRepairs, distinctLoadRepairs, loadRepair } from '../utils/load-repairs.js';
import { t, tm } from '../utils/i18n.js';
import { toRaw } from '../utils/vue-raw.js';
import { specificRulesFileCaptions } from './file-imports.js';
import { isPlainObject, shapeViolations } from '../utils/contract-shape.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import {
  state, sessionOperationAvailability, normalizeLinearSeqList, collapseEmptyLinearSeqList
} from '../state.js';
import { normalizeDefaultColor, normalizeOptionalHexColor, normalizePaletteColors, resolveColorToHex } from '../utils/color-utils.js';
import {
  captureRightDrawerState,
  resetRightDrawerState,
  restoreRightDrawerState
} from './right-drawer-state.js';
import { resetLayoutState, resetSettings as resetSettingsState } from './reset.js';
import { floatSettingsSection, returnSettingsSection } from './settings-pane-state.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { orderedOptionalConservationFiles, parseConservationLabelText, withoutConservationSeries } from './conservation-series.js';
import { cloneJsonData, cloneJsonValue } from './json-clone.js';
import { isNoStrokeColor, normalizeLegendColorOverrides, normalizeStrokeOverrideMap, normalizeStrokeWidth } from './style-overrides.js';
import {
  applyCircularTrackOrderPlacements,
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec
} from './circular-track-slot-model.js';
import {
  applyLinearTrackOrderPlacements,
  LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema
} from './linear-track-slot-model.js';
import { depthFileSlotsFromValue, fitDepthStateToSources, padDepthFileSlots } from './depth-track-state.js';
import { draftDrawsDepth } from './track-legend-rows.js';
import { isUnsetValue, noteReplacedLeaf, truthy } from '../utils/field-reading.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import {
  configFieldRepairs, defaultSettingsContainers, installStoredSettings, readBlacklistText, losatExecutionRepairs, readAnnotationSets, readCanvasPadding,
  readGeneratedCircularPlotTitlePosition, readImportedComparisonAction, draftLinearRecordGap, readLinearRecordGap, readLinearRecordLayoutEnabled,
  readLinearRecordRows, readFeaturePanelTab, readList, readLosatProgram, readPaletteColors, readPlainObject,
  readColorRules, readFileImports, readFilterMode, readQualifierPriorityRules, readRichFeaturePopup, readSessionUiMode, readStringMap, readWhitelistRules,
  featurePlacementsReader, notedLossRepairs, readFeatureOverrides, readPalette, recordsOffReader, resolvedTrackStacks, sessionUiFieldRepairs, sliceFieldRepairs,
  unreadFieldRepairs, valueRepairs, LOSAT_EXECUTION_ROWS, readDownloadDpi, readFeatureColorEdits, readLabelOverrideRows, readLinearSeqFields,
  readSessionInputType, readSimilarityAlignment
} from './config-field-readers.js';
import { drawnLinearSequences, isRecordDrawKey } from './record-draw-selection.js';
import { migrateLegacyLinearLabelVisibility } from './linear-label-visibility.js';
import { migrateLegacyOrthogroupMembers } from './legacy-similarity-alignment.js';
import {
  drawnLegendSideOfMode,
  legendPositionOrDefault,
  readDrawnLegendSideOfMode,
  readLegendPositionOfMode,
  storedLegendPositionRepairs,
  normalizeCircularPlotTitlePosition,
  normalizeLayoutPreferences,
  normalizeLinearPlotTitlePosition,
  replaceLayoutPreferences,
  resolveActiveLayoutPreference,
  restoredLayoutPreferences
} from './layout-preferences.js';
import {
  serializeFeatureVisibilityRules,
  normalizeFeatureVisibilityRule,
  splitLegacyVisibilityRules
} from './feature-visibility.js';
import { canonicalFeatureOverrides, featureDraftMap } from './feature-placement.js';
import {
  ANNOTATION_TARGET_MIGRATION_NOTICE,
  RENDERED_ID_FEATURE_EDIT_FIELDS,
  hasHashColorRules,
  hasRenderedIdFeatureEdits,
  legacyLinearRequestRecords,
  migrateSessionAnnotationTargets,
  migrateSessionFeatureEdits,
  migrateSessionFeaturePlacements,
  migrateSessionHashRules
} from './feature-edit-migration.js';
import {
  buildSessionFeatureRecoveryPlan,
  extractSessionSourceFeatures
} from './session-feature-recovery.js';
import {
  analyzeCatalogSequenceSourceCoverage,
  buildRestoredMatchSequenceSources,
  resolveCircularComparisonSequenceAvailability
} from './match-sequences.js';
import {
  adoptPromotedPreservedOverrides,
  applyLiveSettingWrites,
  buildCanonicalRenderRequest,
  canonicalLinearRecordLayout,
  draftFieldsForConfigValues,
  featureEditsOverProjection,
  linearDefinitionVisibilityOf,
  liveSettingApplies,
  managedConfigOverridePathsForMode,
  projectCanonicalSessionRequest,
  projectSettingsOnlySession,
  readCanonicalResourceRecordCount,
  readStoredRenderRequest
} from './session-request.js';
import {
  createDefaultLinearComparisonPlan,
  linearComparisonEdgeKey,
  losatEdgeOfPair,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan,
  setComparisonEdgeFile
} from './linear-comparisons.js';
import { foldLoadedLinearSources, groupLinearSourceRecords, restoreFoldedLinearSources } from './linear-sources.js';
import { buildSessionResources as assembleSessionResources } from './session-resources.js';
import { ARTIFACT_SLOT_KEYS, createArtifactSlot } from './artifact-slot.js';
import {
  base64ToBytes,
  bytesToBase64,
  getSessionResourceSource,
  readFileBytes
} from './file-content-cache.js';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView,
  isSessionResourceFileView
} from './session-resource-backing.js';
import { normalizeLogicalResults } from './result-normalization.js';
import {
  admitCurrentSessionResults,
  admitLegacyImportedResults,
  createCurrentSessionResultSource,
  createLegacyImportResultSource,
  createSavedResultPlan,
  isCommittedSvgResult
} from './svg-result-ingestion.js';
import {
  admitFeatureCatalog,
  featureStateFromCatalog,
  isAdoptedFeatureCatalog,
  validateFeatureCatalog,
  validateFeatureCatalogForImport
} from './feature-catalog.js';
import { bindRecordDisplayDrafts, migrateLegacyRecordDisplayDrafts, readRecordDisplayDrafts } from './record-display-model.js';
import {
  buildOrthogroupFeatureIndex,
  enrichFeaturesWithOrthogroups,
  normalizeOrthogroupDormantOverrides
} from './orthogroup-feature-metadata.js';
import {
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';
import {
  migrateLegacyLinearComparisonDraft,
  promoteReadGallerySession
} from './gallery-session-migration.js';
import {
  compressSessionData,
  largeSessionDownloadMegabytes
} from './session-file.js';
import { importSessionFile } from './session-import-client.js';
import {
  convertMainSessionComparisonFrames, convertStoredProteinResultFrames, withConvertedResources
} from './main-session-comparison-frame.js';
import { askBeforeLargeDownload, downloadBlob } from './text-download.js';
import { normalizeAnnotationSets } from './annotation-state.js';
import { applySpecificRuleProvenance, savedFeatureColorEdits, savedRuleRowRekeys, savedSpecificRules } from './specific-color-rules.js';
import { legendEntryKey, rekeyLegendRecords } from './legend-svg.js';
import { normalizeLegacyLegendEntryGroups } from './svg-result-normalization.js';
import {
  buildValidatedProteinIdentityIndex,
  classifyRawLosatCacheEntry,
  createLegacyProteinCandidateEnvelope,
  emptyProteinIdentityManifest,
  isCurrentRawLosatCacheEntry,
  isLosatDerivedCacheEntry,
  normalizeLegacyProteinCandidateEnvelope,
  releaseValidatedProteinIdentityIndex,
  serializableLegacyProteinCandidateEnvelope,
  validateDerivedProteinReferences,
  validateProteinRawEntryReferences,
  validateProteinIdentityManifest
} from './losat-cache.js';
import {
  defaultFeatureRendering,
  normalizeFeatureRenderingMap,
  normalizedArrowGeometryState
} from '../utils/feature-rendering.js';
import {
  adoptCurrentSessionDocument,
  adoptRuntimeCanonicalSession,
  isAdoptedCanonicalSession,
  isAuthorityBand,
  isCanonicalBand,
  hasBiologicalSessionInputs,
  holdsNothing,
  hasSettingsOnlyForm,
  isSettingsOnlySessionDocument,
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
  readSessionAuthorityParts,
  readSessionDocument,
  readSessionInput,
  restoresNothing,
  settingsOnlySessionCandidate,
  validateSessionAuthorityInventory
} from './session-authority.js';
import { assertSafeObjectKeysForImport } from './safe-object-keys.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';
import { setResourcePayloadOwner } from './resource-payload-owner.js';
import { migratePersistedWebStateFieldNames, requireCurrentConfigFields } from './current-option-values.js';
import {
  readSimilarityAlignmentResetReceipt,
  validateSimilarityAlignmentResetReceipt,
  CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS,
  LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
  reconcileImportedLinearTypographyLink,
  holdsCliWriterConfig,
  requireActiveConfigContainers,
  validateCurrentWriterActiveConfig,
  draftFeatureOverlapTolerance
} from './session-active-config-contract.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS,
  classifyImportedComparisonIntent,
  createImportedComparisonIntentState,
  restoreImportedComparisonIntent,
  serializeImportedComparisonResolution
} from './imported-comparison-intent.js';

import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import {
  CANONICAL_SESSION_MIN_VERSION,
  CURRENT_ARTIFACT_SESSION_MIN_VERSION,
  CURRENT_AUTHORITY_SESSION_MIN_VERSION,
  DRAWN_HASH_RULE_SESSION_MAX_VERSION,
  FILE_IMPORT_ORIGIN_SESSION_VERSION,
  LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
  LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  LOSAT_DERIVED_CACHE_SCHEMA,
  MODE_SCOPED_SESSION_VERSION,
  ROW_KEYED_LEGEND_SESSION_VERSION,
  SESSION_VERSION,
  TYPED_DRAFT_SESSION_VERSION
} from '../contracts/schema.generated.js';
import {
  SLICE_MODES,
  sessionDepthSourceWidths,
  sessionLegendEditsByRowKey,
  splitDraftIntoModes,
  unscopedDraftRows,
  validateModeSliceFields
} from './mode-scoped-migration.js';

const { nextTick, ref } = window.Vue;

// How LOSAT runs: one app-level setting for both drawings (Session 46
// `ui.losatExecution`; a Session 44 or older draft keeps it in `config.losat`).
const LOSAT_EXECUTION_FIELDS = MODE_SCOPED_SETTINGS.losatExecutionFields;
/**
 * Installs saved LOSAT execution settings; a field the value lacks takes its
 * normalized default, as a saved draft always read.
 * @param {SessionLosatExecution} source Unvalidated: each field is matched or converted here.
 */
export const applyLosatExecutionData = (source) => {
  LOSAT_EXECUTION_ROWS.forEach(({ field, read }) => {
    state.losatExecution[field] = read(source[/** @type {keyof SessionLosatExecution} */ (field)], source).value;
  });
};
/** @returns {Record<string, any>} The settings as Save writes them (`ui.losatExecution`). */
export const buildLosatExecutionData = () => cloneJsonData(Object.fromEntries(
  LOSAT_EXECUTION_FIELDS.map((/** @type {string} */ field) => [field, state.losatExecution[field]])
));

/** @import { ActiveWebConfig, AlignmentResetCanonical } from './session-active-config-contract.js' */
/** @import { SettingLoss, SettingsContainers } from './config-field-readers.js' */
/** @import { LegendEditsUnread, StrokeOverrideLosses } from './style-overrides.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { StoredDrawing } from './main-session-comparison-frame.js' */
/** @import { CommittedCanonicalSession, LiveSettingWrite, RequestCanonicalComparison, SessionRequestProjection,
 *   SettingsOnlyProjectionInput } from './session-request.js' */
/** @import { ImportedComparisonIntent } from './imported-comparison-intent.js' */
/** @import { FeatureCatalog, FeatureCatalogResult } from './feature-catalog.js' */
/** @import { LegacyResultNormalization, SavedResultEdits, SvgResultEntry } from './svg-result-ingestion.js' */
/** @import { DroppedLegendRecord, PythonLegendKey, PythonLegendRow } from './legend-svg.js' */
/** @import { SessionResourceSource } from './session-resources.js' */
/** @import { LoadedTrackSettings } from './track-legend-rows.js' */
/** @import { ArtifactSlot } from './artifact-slot.js' */
/** @import { LinearComparisonEdge, LinearComparisonPlan, LinearComparisonResolvedEdge } from './linear-comparisons.js' */
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy } from '../state.js' */
/** @import { LinearRecordRow } from './linear-record-layout.js' */
/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { StoredLayoutPreferences } from './layout-preferences.js' */
/** @import { TrackSlotRow } from './track-slot-validation.js' */
/** @import { DiagnosticContext, UserFacingError } from '../utils/error-normalization.js' */
/** @import { ConfirmLoadRepairs, LoadRepair, Repaired } from '../utils/load-repairs.js' */
/** @import { ProteinIdentityManifest } from './losat-cache.js' */
/** @import { DownloadChoicePort } from './text-download.js' */

// The Session document's types are generated from gbdraw/session_contract.py:
// `GbdrawSession` is the writer's type (the current version only); a reader
// takes unvalidated data (R14).
/** @import { GbdrawSession, SessionChecked, SessionChecked27, SessionChecked31, SessionChecked39, SessionEnvelope, SessionInput, SessionLosatExecution, SessionModeSlice, SessionOtherModeResult } from '../contracts/session.generated.js' */
/** @import { ResultSetPath, StoredRequestSchemas, UnvalidatedFields } from './session-authority.js' */
/** @import { FileLike, RawSessionResourceDescriptor, SessionResourceFileView } from './session-resource-backing.js' */
/** @import { DepthFilePayload } from './depth-file-codec.js' */
/** @import { RawConservationSeriesEntry } from './conservation-series.js' */
/** @import { DepthFileSlots } from './depth-track-state.js' */

/** @import { FallbackUiData, FeatureEditsData, OffsetData } from './history-snapshot.js' */

// A saved `ui` or `features` object (a Session's, a mode slice's, a History
// snapshot's) has one unvalidated view, declared with History's fallback
// readers (services/history-snapshot.js); these name the parts read here.
/**
 * The saved palette values of a `ui` object as they are read: unvalidated.
 * @typedef {Pick<FallbackUiData, 'appliedPaletteName' | 'appliedPaletteColors' | 'pendingPaletteName'
 *   | 'pendingPaletteColors'>} SavedPaletteFields
 */
/**
 * A saved `ui` object as `applyUiStateData` reads it: every leaf is
 * unvalidated until the line that reads it converts or matches it.
 * @typedef {Omit<FallbackUiData, 'selectedResultIndex' | 'layoutPreferences'> & {
 *   losatExecution?: unknown, featurePanelTab?: unknown,
 *   circularPlotTitlePosition?: unknown, zoom?: unknown, canvasPan?: OffsetData | null
 * }} SavedUiFields
 */
/**
 * Saved feature state (a Session's `features`, a mode slice's, a History
 * checkpoint's) as it is read: unvalidated.
 * @typedef {FeatureEditsData & { extractedFeatures?: unknown, biologicalFeatures?: unknown,
 *   featureRecordIds?: unknown, labelOverrideRows?: unknown }} SavedFeatureFields
 */

/**
 * What a loaded Result needs from the Session: its mode and, for the Result
 * set an older Session's legacy composition applies to, the committed legend
 * and title sides with the saved user offsets (null for another set). Its
 * saved bytes already show its strokes.
 * @typedef {object} LegacyResultSvgData
 * @property {{ legendSide: string, titleSide: string, userDeltas: Record<string, number[] | null> } | null} composition
 * @property {'circular' | 'linear'} mode
 * @property {LoadRepair[] | null} repairs Collects what the transform cannot restore, for the Load dialog.
 * @property {string} path Where the Result's set is in the file.
 */

/**
 * The composition root's transform of a loaded Result (R13 port): it gives a
 * Result without composition metadata the legacy composition, and a Result
 * with metadata a legend side of its mode (OV-672). Returns whether the SVG
 * changed. `appliesToContent` tells from a Result's text whether it lacks
 * composition metadata (and `composes`, its set has a legacy composition) or
 * stores a side outside `mode`'s positions, so a Session 40+ Result is parsed
 * only when it does (OV-273).
 * @typedef {((svg: SVGSVGElement, data: LegacyResultSvgData) => boolean)
 *   & { appliesToContent: (content: unknown, mode: 'circular' | 'linear', composes: boolean) => boolean }} LegacyResultSvgTransform
 */
export { SESSION_VERSION };
const LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32;
const LOSAT_DERIVED_CACHE_LIMIT = 16;
// D-25 (PD-OI-079): a Result without current feature metadata (a legacy
// Session) is saved only after one Generate; the error offers that Generate.
// E1: `diagramMode` names the mode whose Result it is when that mode is not
// shown; the error's Generate action then runs in that mode.
/** @param {Record<string, any>} [context] */
const sessionSaveRequiresGenerate = (context = {}) => diagnosticError(
  'SESSION_SAVE_REQUIRES_GENERATE', context, /** @type {Record<string, any>} */ ({
    operation: 'session-save', stage: 'result-admission'
  })
);
// A Save wrapper keeps a recognized cause; otherwise it reports the bounded fallback.
/**
 * @template F
 * @param {Error | string | Record<string, unknown> | null | undefined} error
 * @param {F} fallback
 */
const recognizedCauseOr = (error, fallback) => {
  const model = normalizeUserFacingError(error);
  return model && !['UNKNOWN', 'VALIDATION_UNCLASSIFIED'].includes(model.code)
    ? diagnosticError(model.code, model.context, { stage: model.stage }) : fallback;
};

/**
 * @template {UnvalidatedFields<string>} T
 * @param {T} colors A color map. The checker gives the copy the type `{}` for every caller, so the copy does not
 *   carry the type of `colors`: a copy of a saved, unvalidated map is still unvalidated.
 */
const cloneColors = (colors) => ({ ...(colors || {}) });

/** @param {unknown} colors */
const hasColorEntries = (colors) =>
  Boolean(colors) && typeof colors === 'object' && !Array.isArray(colors)
    // The cast: `Boolean(colors)` above proves the object is not null.
    && Object.keys(/** @type {NonNullable<UnvalidatedFields<string>>} */ (colors)).length > 0;

/** @param {unknown} colors */
const normalizeColorMap = (colors) => {
  /** @type {Record<string, string>} */
  const normalized = {};
  if (!colors || typeof colors !== 'object' || Array.isArray(colors)) return normalized;
  Object.entries(colors).forEach(([key, value]) => {
    normalized[key] = resolveColorToHex(String(value || '').trim());
  });
  return normalized;
};

/** @param {unknown} paletteName */
const paletteColorsFromDefinitions = (paletteName) => {
  const name = String(paletteName || '').trim();
  if (!name) return null;
  /** @type {Record<string, Record<string, string>>} */
  const definitions = state.paletteDefinitions?.value || {};
  const colors = definitions[name];
  return hasColorEntries(colors) ? normalizePaletteColors(cloneColors(colors)) : null;
};

/** @param {unknown} source */
const cloneStringMap = (source) => readStringMap(source).value;

/** @param {unknown} rules */
const normalizeFeatureVisibilityRulesForSession = (rules) => (
  Array.isArray(rules) ? rules.map((rule) => normalizeFeatureVisibilityRule(rule)) : []
);

// A Session before manual rules had their own field stored the editor's
// per-feature rows among the rules (`featureVisibilityRules`); its rows keyed
// by rendered ID join the other rendered-ID edits for migration.
/** @param {Record<string, unknown>} [features] A stored feature-edit draft. @returns {Record<string, unknown>} */
const splitLegacyFeatureVisibilityRules = (features = {}) => {
  if (Array.isArray(features.featureVisibilityManualRules) || !Array.isArray(features.featureVisibilityRules)) {
    return features;
  }
  const { manualRules, overrides } = splitLegacyVisibilityRules(features.featureVisibilityRules);
  return {
    ...features,
    featureVisibilityManualRules: manualRules,
    featureVisibilityOverrides: { ...overrides, ...(features.featureVisibilityOverrides || {}) }
  };
};

// The identity-keyed per-feature edit draft: the rows a request accepts
// (services/feature-placement.js); Load lists the others (`sliceFieldRepairs`).
/** @param {unknown} value */
const featureOverridesForState = (value) => featureDraftMap(
  canonicalFeatureOverrides(readFeatureOverrides(value).value)
);

/**
 * @param {UnvalidatedFields<'nucleotide_sequence' | 'amino_acid_sequence' | 'nucleotideSequence'
 *   | 'aminoAcidSequence'>} feature
 */
const sanitizeExtractedFeatureForSession = (feature) => {
  if (!feature || typeof feature !== 'object' || Array.isArray(feature)) return feature;
  const {
    nucleotide_sequence: _nucleotideSequence,
    amino_acid_sequence: _aminoAcidSequence,
    nucleotideSequence: _nucleotideSequenceAlias,
    aminoAcidSequence: _aminoAcidSequenceAlias,
    ...rest
  } = feature;
  return rest;
};

/** @param {unknown} features */
const sanitizeExtractedFeaturesForSession = (features) => {
  if (!Array.isArray(features)) return [];
  return features.map((feature) => sanitizeExtractedFeatureForSession(feature));
};

/**
 * @param {Record<string, string>} target
 * @param {unknown} source
 */
const replaceStringMap = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(cloneStringMap(source)).forEach(([key, value]) => {
    target[key] = value;
  });
};


/**
 * @param {DrawingState} drawing
 * @param {unknown} rules
 */
const replaceQualifierPriorityRules = (drawing, rules) => {
  drawing.manualPriorityRules.splice(
    0,
    drawing.manualPriorityRules.length,
    ...readQualifierPriorityRules(rules).value
  );
};

/** @param {unknown} value */
const parseMultiRecordPositionToken = (value) => {
  const raw = String(value ?? '').trim();
  const separatorIndex = raw.lastIndexOf('@');
  if (separatorIndex <= 0 || separatorIndex === raw.length - 1) return null;

  const selector = raw.slice(0, separatorIndex).trim();
  const row = Number(raw.slice(separatorIndex + 1).trim());
  if (!selector || !Number.isInteger(row) || row <= 0) return null;
  return { selector, row };
};

/**
 * One row of `adv.multi_record_positions`: a record selector and its 1-based row on the Multi-Record Canvas.
 * `applyConfigData` writes the checked rows.
 * @typedef {{ selector: string, row: number }} MultiRecordPosition
 */
/** @param {{ args?: unknown } | null | undefined} cliInvocation */
const multiRecordPositionsFromCliInvocation = (cliInvocation) => {
  const args = Array.isArray(cliInvocation?.args) ? cliInvocation.args : [];
  /** @type {MultiRecordPosition[]} */
  const positions = [];
  const seenSelectors = new Set();

  args.forEach((arg, index) => {
    if (arg !== '--multi_record_position') return;
    const position = parseMultiRecordPositionToken(args[index + 1]);
    if (!position || seenSelectors.has(position.selector)) return;
    seenSelectors.add(position.selector);
    positions.push(position);
  });

  return positions;
};

/**
 * @param {Record<string, unknown>} config A restored draft; its `adv` is written.
 * @param {{ args?: unknown } | null | undefined} cliInvocation
 */
const hydrateMissingMultiRecordPositionsFromCliInvocation = (config, cliInvocation) => {
  if (!isPlainObject(config) || !isPlainObject(config.form) || config.form.multi_record_canvas !== true) {
    return;
  }

  const adv = isPlainObject(config.adv) ? config.adv : {};
  if (Array.isArray(adv.multi_record_positions)) return;

  const positions = multiRecordPositionsFromCliInvocation(cliInvocation);
  if (positions.length === 0) return;

  adv.multi_record_positions = positions;
  config.adv = adv;
};

/**
 * @param {Record<string, any> | null | undefined} invocation
 * @returns {boolean}
 */
const isCliInvocationSessionExportable = (invocation) => {
  if (!invocation || typeof invocation !== 'object') return false;
  if (invocation.sessionExportable === false) return false;
  const bindings = Array.isArray(invocation.fileBindings) ? invocation.fileBindings : [];
  return bindings.every((binding) => String(binding?.slot || '').startsWith('files.'));
};

/** @param {unknown} name */
const makeSafeFilename = (name) => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || 'gbdraw_session';
};

/** @param {unknown} title */
const buildSessionFilename = (title) => {
  const base = String(title || '').trim();
  if (!base) return 'gbdraw_session.json.gz';
  const safe = makeSafeFilename(base);
  return `${safe}.gbdraw-session.json.gz`;
};

// Load reads a stored configuration's older Linear visibility (every caller is a Load path).
/** @template T @param {T} configData @returns {T} */
const withCurrentLinearLabelVisibility = (configData) => {
  if (!isPlainObject(configData) || !isPlainObject(configData.adv)) return configData;
  // The cast: the same draft with its `adv` read again.
  return /** @type {T} */ ({
    ...configData,
    adv: migrateLegacyLinearLabelVisibility(configData.adv, { load: true })
  });
};

/**
 * @param {Record<string, any>} [configData] The stored active configuration.
 * @returns {Record<string, any>}
 */
const migrateImportedCircularTrackSlots = (configData = {}) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return configData;
  // Only schema-3 slots migrate: the stack's reader reads any other stack in
  // the current schema and lists what it cannot use (OV-758).
  if (!Array.isArray(adv.circular_track_slots)
    || adv.circular_track_slots_schema_version !== LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION) return configData;

  const defaultNt = adv.nt || 'GC';
  const preset = configData.form?.track_type || 'tuckin';
  return {
    ...configData,
    adv: {
      ...adv,
      circular_track_slots_schema_version: CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
      circular_track_slots: adv.circular_track_slots.map((/** @type {string | TrackSlotRow} */ slot, /** @type {number} */ index) => (
        typeof slot === 'string'
          ? parseCircularTrackSlotSpec(
              migrateLegacyCircularTrackSlotSpec(slot),
              index,
              defaultNt,
              preset
            )
          : migrateLegacyCircularTrackSlot(slot)
      ))
    }
  };
};

const migratePersistedWebOptionValues = (configData = {}) => {
  if (!configData || typeof configData !== 'object' || Array.isArray(configData)) {
    return configData;
  }
  // The older option values (`spreadout`, `on_feature`, `sqrt`) are read as
  // their current names by the reader table, at every version.
  return withCurrentLinearLabelVisibility(migratePersistedWebStateFieldNames(configData));
};

/**
 * @param {Record<string, any>} [configData] The stored active configuration.
 * @param {number | null} [sourceSessionVersion]
 * @returns {Record<string, any>}
 */
const migrateImportedLinearTrackSlots = (configData = {}, sourceSessionVersion = null) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return configData;
  // A stack of another kind or schema is not migrated: its reader reads it in
  // the current schema and lists what it cannot use (OV-758).
  if (!Array.isArray(adv.linear_track_slots)) return configData;
  const storedSchemaVersion = adv.linear_track_slots_schema_version;
  const sessionUsesLegacySemantics = (
    Number.isInteger(sourceSessionVersion) &&
    // Number.isInteger above is false for null, so the version is a number here.
    /** @type {number} */ (sourceSessionVersion) <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
  );
  const sourceSchemaVersion = sessionUsesLegacySemantics
    ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
    : storedSchemaVersion;
  if (![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(sourceSchemaVersion)) return configData;

  return {
    ...configData,
    adv: {
      ...adv,
      linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
      linear_track_slots: migrateLinearTrackSlotsToCurrentSchema(
        adv.linear_track_slots,
        sourceSchemaVersion
      )
    }
  };
};

/** @param {unknown} value */
const hasStoredLayoutValue = (value) => typeof value === 'string' && value.trim() !== '';

// A stored entry color is read with the Default colors domain (D-41). An entry
// whose color lies outside it is dropped and, when `repairs` is given, listed
// there for the Load dialog (D-43) under `path`, the stored Legend. A list
// holds one entry per row key (Session 47: two rows may show one name); a
// Session 46 or older kept the first entry of each caption (OV-719). A later
// entry of one row is dropped and listed, and so is an entry that is no object
// or has no caption (TB8).
/**
 * @param {unknown} entries
 * @param {'entries' | 'dormantEntries' | 'deletedEntries'} list
 * @param {{ repairs: LoadRepair[] | null, path: string, sessionVersion: number }} load
 */
const normalizeSessionLegendEntries = (entries, list, { repairs, path, sessionVersion }) => {
  /** @type {SettingLoss[]} */
  const losses = [];
  if (!Array.isArray(entries)) {
    if (!isUnsetValue(entries)) repairs?.push(...notedLossRepairs([{ at: list, stored: entries, used: [] }], path));
    return [];
  }
  const rowOf = sessionVersion >= ROW_KEYED_LEGEND_SESSION_VERSION
    ? (/** @type {Record<string, unknown>} */ entry) => legendEntryKey(entry)
    : (/** @type {Record<string, unknown>} */ entry) => String(entry.caption || '').trim();
  /** @type {Record<string, unknown>[]} */
  const normalized = [];
  /** @type {Set<string>} */
  const keys = new Set();
  /** @type {{ list: string, caption: string, key: string }[]} */
  const dropped = [];
  /** @type {Map<string, number>} An entry that holds something but is no object or has no caption, by list. */
  const unread = new Map();
  /** @param {string} from */
  const countUnread = (from) => unread.set(from, (unread.get(from) ?? 0) + 1);
  entries.forEach((entry, index) => {
    if (!isPlainObject(entry)) {
      if (!isUnsetValue(entry)) countUnread(list);
      return;
    }
    const caption = String(entry.caption || '').trim();
    const key = rowOf(entry);
    const color = normalizeDefaultColor(entry.color);
    const from = entry.dormant === true ? 'dormantEntries' : list;
    if (!caption && Object.values(entry).some((field) => !isUnsetValue(field))) countUnread(from);
    if (caption && !color && !keys.has(key)) dropped.push({ list: from, caption, key });
    if (caption && color && keys.has(key)) {
      repairs?.push(loadRepair({ code: 'legend-entry-duplicate', path: `${path}.${from}`, effect: 'dropped', label: caption }));
    }
    if (!caption || !color || keys.has(key)) return;
    keys.add(key);
    // `showStroke`, the Stroke options disclosure that earlier Sessions saved,
    // is view state (OV-157).
    const stored = cloneJsonData(entry);
    delete stored.showStroke;
    if (!Array.isArray(entry.featureIds) && !isUnsetValue(entry.featureIds)) {
      losses.push({ at: `${list}[${index}].featureIds`, stored: entry.featureIds, used: [] });
    }
    normalized.push({
      ...stored,
      caption,
      color,
      featureIds: Array.isArray(entry.featureIds)
        ? entry.featureIds.map((/** @type {unknown} */ id) => String(id || '').trim()).filter(Boolean)
        : []
    });
  });
  repairs?.push(...notedLossRepairs(losses, path));
  unread.forEach((count, from) => repairs?.push(loadRepair({
    code: 'setting-entries-dropped', path: `${path}.${from}`, effect: 'dropped', count
  })));
  // A row that a later entry with a valid color keeps is shown, not dropped.
  dropped.filter(({ key }) => !keys.has(key)).forEach(({ list: from, caption }) => repairs?.push(loadRepair({
    code: 'legend-color-dropped', path: `${path}.${from}`, effect: 'dropped', label: caption
  })));
  return normalized;
};

/**
 * @template T
 * @param {T[] | null | undefined} value An app-written list (a state value).
 * @returns {T[]}
 */
const cloneJsonArray = (value) => {
  if (!Array.isArray(value)) return [];
  return cloneJsonValue(value, []);
};

/**
 * @param {unknown} value
 * @returns {Record<string, unknown>} A copy of a plain object, or an empty one.
 */
const cloneJsonObject = (value) => {
  if (!isPlainObject(value)) return {};
  return cloneJsonValue(value, {});
};

/**
 * @param {unknown} source
 * @param {string} at Where `source` is, for `losses`.
 * @param {SettingLoss[] | null} [losses] Written: a list that is no list, an entry read as other text, and the
 *   number of entries that hold something but read as no text.
 */
const normalizeStringArray = (source, at, losses = null) => {
  if (!Array.isArray(source)) {
    if (!isUnsetValue(source)) losses?.push({ at, stored: source, used: [] });
    return [];
  }
  let dropped = 0;
  const read = source.flatMap((value, index) => {
    const text = String(value || '').trim();
    if (!text && !isUnsetValue(value)) dropped += 1;
    else if (text) noteReplacedLeaf(losses, `${at}[${index}]`, value, text);
    return text ? [text] : [];
  });
  if (dropped) losses?.push({ at, stored: source, used: read, dropped });
  return read;
};

/** @type {string | null} */
let lastSessionFilename = null;
/** @type {Record<string, any> | null} */
let preservedCliOptions = null;
// The variable keeps its earlier type: `applyResultsData` passes the request's
// `diagramOptions.selectedFeaturesSet` on as feature types without a check,
// which the true type refuses (stage 1d).
/** @type {Record<string, any> | null} */
let committedCanonicalSession = null;
// Counts replacements of the committed Session, so a computed that reads the
// committed request follows Generate, Load, a History restore, and a mode switch.
const committedCanonicalRevision = ref(0);
/** @param {Record<string, any> | null} next */
const setCommittedCanonicalSession = (next) => {
  committedCanonicalSession = next;
  committedCanonicalRevision.value += 1;
};
// OV-484: the track settings a Session without a request (a Session 27 to 30)
// was saved with, which stand in for the request that drew its Result (they are
// the saved draft). Load is the one writer, and only a Load replaces or
// clears the record; it is read while the Result's mode has no committed
// request, so Generate, a History step and a mode switch leave it as it is.
/** @type {LoadedTrackSettings | null} */
let loadedTrackSettings = null;
/** @type {ReturnType<typeof adoptCurrentSessionResources> | null} */
let activeSessionResourceTable = null;
/** @type {Record<string, any> | null} */
let adoptedProteinIdentityManifest = null;
// Current-session preflight receipts pair each adopted cache value with the
// exact manifest that already validated its protein references and raw text.
const adoptedLosatCacheValues = new WeakMap();

/** @param {{ renderRequest?: unknown, resources?: unknown, webFiles?: unknown } | null | undefined} canonical */
const cloneCanonicalSession = (canonical) => {
  if (
    !canonical
    || !isPlainObject(canonical.renderRequest)
    || !isPlainObject(canonical.resources)
  ) return null;
  return {
    renderRequest: cloneJsonData(canonical.renderRequest),
    resources: cloneJsonData(canonical.resources),
    webFiles: isPlainObject(canonical.webFiles)
      ? cloneJsonData(canonical.webFiles)
      : {}
  };
};

/**
 * @param {LinearComparisonPlan} target
 * @param {unknown} source
 */
const replaceLinearComparisonPlan = (target, source) => {
  const normalized = normalizeLinearComparisonPlan(source);
  target.mode = normalized.mode;
  target.defaultSource = normalized.defaultSource;
  if (!Array.isArray(target.edges)) target.edges = [];
  target.edges.splice(0, target.edges.length, ...normalized.edges);
};

/** @param {unknown} plan */
const serializeLinearComparisonPlan = (plan) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  return {
    mode: normalized.mode,
    defaultSource: normalized.defaultSource,
    edges: normalized.edges.map((edge) => ({
      id: edge.id,
      queryUid: edge.queryUid,
      subjectUid: edge.subjectUid,
      included: edge.included,
      fileActive: edge.fileActive,
      losatFilenameActive: edge.losatFilenameActive,
      source: edge.source,
      losatFilename: edge.losatFilename
    }))
  };
};

// A drawing's settings as its Session 46 slice `config` holds them: the
// registry rows of `config` (mode-scoped-settings.generated.js) and nothing
// else. App-level settings (LOSAT execution, Instant Preview, the rich popup)
// and the Session's CLI provenance are saved outside the slices.
/** @param {DrawingState} drawing */
export const buildConfigData = (drawing) => ({
  form: drawing.form,
  adv: {
    ...drawing.adv,
    feature_overlap_tolerance_bp: draftFeatureOverlapTolerance(drawing.adv),
    feature_shapes: {
      repeat_region: defaultFeatureRendering('repeat_region'),
      ...normalizeFeatureRenderingMap(drawing.adv.feature_shapes || {})
    }
  },
  losat: cloneJsonData(drawing.losat || {}),
  colors: drawing.currentColors.value,
  palette: drawing.selectedPalette.value,
  rules: drawing.manualSpecificRules,
  qualifierPriorityRules: readQualifierPriorityRules(drawing.manualPriorityRules).value,
  filterMode: drawing.filterMode.value,
  whitelist: drawing.manualWhitelist,
  blacklistText: drawing.manualBlacklist.value,
  fileImports: cloneJsonData(drawing.fileImports),
  losatProgram: drawing.losatProgram.value,
  circularConservation: drawing.circularConservation,
  annotationSets: normalizeAnnotationSets(drawing.annotationSets),
  recordDisplayDrafts: cloneJsonData(drawing.recordDisplayDrafts),
  // Written only when a record is OFF: a Session with every record ON keeps its bytes (omission = all ON).
  ...(drawing.recordsOff?.length ? { recordsOff: [...drawing.recordsOff] } : {}),
  featurePlacementOverrides: cloneJsonData(drawing.featurePlacementOverrides),
  linearRecordLayout: {
    enabled: Boolean(drawing.linearRecordLayoutEnabled.value),
    recordGap: Number(draftLinearRecordGap(drawing.linearRecordGap.value)) || 0,
    // The rows as Load reads them (one reader; a row state holds is valid).
    rows: readLinearRecordRows(drawing.linearRecordRows).value
  },
  linearComparisonPlan: serializeLinearComparisonPlan(drawing.linearComparisonPlan),
  importedComparisonResolution: serializeImportedComparisonResolution(
    drawing.importedComparisonIntent
  ),
  unmanagedConfigOverrides: cloneJsonData(drawing.unmanagedConfigOverrides || {}),
  webEdits: {
    orthogroupNameOverrides: cloneStringMap(drawing.orthogroupNameOverrides),
    orthogroupDescriptionOverrides: cloneStringMap(drawing.orthogroupDescriptionOverrides),
    orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(drawing.orthogroupDormantOverrides)
  }
});

const defaultEditorStateData = () => ({
  legend: {
    entries: [],
    deletedEntries: [],
    dormantEntries: [],
    originalOrder: [],
    originalColors: {},
    colorOverrides: {},
    strokeOverrides: {},
    addedCaptions: []
  },
  featureStrokes: {
    overrides: {}
  },
  originalSvgStroke: {
    color: null,
    width: null
  },
  featureCatalog: null
});

// State, History, and Session rollback hold the Generate-owned catalog by
// reference; only an admitted catalog may enter state. Vue proxies obscure the
// identity that admission records.
/**
 * @param {unknown} catalog A catalog in state, or one a saved draft or snapshot holds.
 * @returns {FeatureCatalog | null}
 */
const admittedFeatureCatalog = (catalog) => {
  const rawCatalog = toRaw(catalog);
  return isAdoptedFeatureCatalog(rawCatalog) ? rawCatalog : null;
};

/** @param {DrawingState} drawing */
export const buildEditorStateData = (drawing) => ({
  legend: {
    entries: cloneJsonArray(drawing.legendEntries.value),
    deletedEntries: cloneJsonArray(drawing.deletedLegendEntries.value),
    dormantEntries: cloneJsonArray(drawing.dormantLegendEntries.value),
    originalOrder: cloneJsonArray(state.originalLegendOrder.value),
    originalColors: cloneStringMap(state.originalLegendColors.value),
    colorOverrides: cloneJsonObject(drawing.legendColorOverrides),
    strokeOverrides: cloneJsonObject(drawing.legendStrokeOverrides),
    addedCaptions: Array.from(drawing.addedLegendCaptions.value || [])
      .map((caption) => String(caption || '').trim())
      .filter(Boolean)
  },
  featureStrokes: {
    overrides: cloneJsonObject(drawing.featureStrokeOverrides)
  },
  originalSvgStroke: {
    color: state.originalSvgStroke.value?.color ?? null,
    width: state.originalSvgStroke.value?.width ?? null
  },
  alignmentResetReceipt: cloneJsonValue(state.similarityAlignmentResetReceipt?.value, null),
  featureCatalog: admittedFeatureCatalog(state.featureCatalog?.value)
});

/**
 * The losses of a stroke override map read at `at` (a dropped entry is
 * counted, and a count names no value).
 * @param {StrokeOverrideLosses} losses
 * @param {string} at
 * @returns {SettingLoss[]}
 */
const strokeOverrideLosses = ({ leaves, dropped }, at) => [
  ...leaves.map((leaf) => ({ ...leaf, at: `${at}.${leaf.at}` })),
  ...(dropped ? [{ at, stored: null, used: {}, dropped }] : [])
];
/**
 * The rows of the Legend edits Load could not read (TB9) of the Legend stored
 * at `path`: the color edits dropped (one row per caption) and the stroke
 * edits' losses.
 * @param {LegendEditsUnread} unread
 * @param {string} path
 * @returns {LoadRepair[]}
 */
const legendEditRepairs = ({ colors, strokes }, path) => [
  ...notedLossRepairs(strokeOverrideLosses(strokes, 'strokeOverrides'), path),
  ...[...new Set(colors)].map((caption) => loadRepair({
    code: 'legend-color-dropped', path: `${path}.colorOverrides`, effect: 'dropped', label: caption
  }))
];

/**
 * @param {unknown} [editorState] A stored editor state; one that is no object is read as none.
 * @param {{ featureCatalog?: unknown, repairs?: LoadRepair[] | null, path?: string, sessionVersion?: number }} [options]
 *   `repairs` collects what Load cannot restore of the editor state stored at `path`; `sessionVersion`: the
 *   version of the Session it was read from.
 */
export const normalizeEditorStateData = (editorState = {}, {
  featureCatalog = undefined, repairs = null, path = 'editorState', sessionVersion = SESSION_VERSION
} = {}) => {
  const defaults = defaultEditorStateData();
  const source = isPlainObject(editorState) ? editorState : {};
  const legend = isPlainObject(source.legend) ? source.legend : {};
  const featureStrokes = isPlainObject(source.featureStrokes) ? source.featureStrokes : {};
  const originalSvgStroke = isPlainObject(source.originalSvgStroke) ? source.originalSvgStroke : {};
  // A Session 46 slice lists the renamed rows its Result does not draw
  // (OV-120) after the shown rows, marked `dormant`.
  const legendLoad = { repairs, path: `${path}.legend`, sessionVersion };
  const entries = normalizeSessionLegendEntries(legend.entries, 'entries', legendLoad);
  const dormantEntries = [
    ...normalizeSessionLegendEntries(legend.dormantEntries, 'dormantEntries', legendLoad),
    ...entries.filter((entry) => entry.dormant === true)
  ].map(({ dormant: _dormant, ...entry }) => entry);
  // What Load cannot read of the edits is listed when `repairs` is given (TB9);
  // an edit that holds nothing (an empty override) is no row.
  /** @type {SettingLoss[] | null} */
  const losses = repairs ? [] : null;
  /** @type {LegendEditsUnread | null} The Legend's color edits and the stroke edits of both stores. */
  const unreadEdits = repairs ? { colors: [], strokes: { leaves: [], dropped: 0 } } : null;
  const droppedColors = unreadEdits?.colors ?? null;
  /**
   * A map of edits that is no object is read as none; an empty value ('',
   * null, [] or {}) holds no edit, so it loses nothing (E-OD36).
   * @param {unknown} map
   * @param {string} at
   */
  const noteNoMap = (map, at) => {
    if (!isUnsetValue(map) && !holdsNothing(map) && !isPlainObject(map)) losses?.push({ at, stored: map, used: {} });
  };
  // A Legend of another kind loses the Result's original Legend order and colors (R2-M5).
  noteNoMap(source.legend, 'legend');
  /**
   * @param {unknown} map
   * @param {string} at
   */
  const readColors = (map, at) => {
    noteNoMap(map, at);
    return normalizeLegendColorOverrides(map, droppedColors);
  };
  /**
   * @param {unknown} map
   * @param {string} at
   */
  const readStrokeOverrides = (map, at) => {
    noteNoMap(map, at);
    /** @type {StrokeOverrideLosses | null} */
    const strokeLosses = losses ? { leaves: [], dropped: 0 } : null;
    const read = normalizeStrokeOverrideMap(map, { requireOverride: true, losses: strokeLosses });
    if (losses && strokeLosses) losses.push(...strokeOverrideLosses(strokeLosses, at));
    return read;
  };
  /**
   * @template T
   * @param {unknown} stored
   * @param {string} field
   * @param {(value: unknown) => T | null} read
   * @param {T | null} fallback
   * @returns {T | null}
   */
  const readSvgStroke = (stored, field, read, fallback) => {
    if (!Object.hasOwn(originalSvgStroke, field)) return fallback;
    const value = read(stored);
    if (!isUnsetValue(stored) && value === null && !isNoStrokeColor(stored)) losses?.push({ at: `originalSvgStroke.${field}`, stored, used: null });
    return value;
  };
  if (!isPlainObject(source.originalSvgStroke) && !isUnsetValue(source.originalSvgStroke)) {
    losses?.push({ at: 'originalSvgStroke', stored: source.originalSvgStroke, used: defaults.originalSvgStroke });
  }

  const read = {
    legend: {
      entries: entries.filter((entry) => entry.dormant !== true),
      deletedEntries: normalizeSessionLegendEntries(legend.deletedEntries, 'deletedEntries', legendLoad),
      dormantEntries,
      originalOrder: normalizeStringArray(legend.originalOrder, 'legend.originalOrder', losses),
      originalColors: readColors(legend.originalColors, 'legend.originalColors'),
      colorOverrides: readColors(legend.colorOverrides, 'legend.colorOverrides'),
      strokeOverrides: readStrokeOverrides(legend.strokeOverrides, 'legend.strokeOverrides'),
      addedCaptions: normalizeStringArray(legend.addedCaptions, 'legend.addedCaptions', losses)
    },
    featureStrokes: {
      overrides: readStrokeOverrides(featureStrokes.overrides, 'featureStrokes.overrides')
    },
    originalSvgStroke: {
      color: readSvgStroke(originalSvgStroke.color, 'color', normalizeOptionalHexColor, defaults.originalSvgStroke.color),
      width: readSvgStroke(originalSvgStroke.width, 'width', normalizeStrokeWidth, defaults.originalSvgStroke.width)
    },
    alignmentResetReceipt: cloneJsonValue(source.alignmentResetReceipt, null),
    featureCatalog: admittedFeatureCatalog(
      featureCatalog !== undefined ? featureCatalog : source.featureCatalog
    )
  };
  if (repairs && losses && unreadEdits) repairs.push(...notedLossRepairs(losses, path), ...legendEditRepairs(unreadEdits, `${path}.legend`));
  return read;
};

/**
 * @param {Record<string, unknown>} target
 * @param {Record<string, unknown> | null | undefined} source
 */
const replacePlainObject = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

/**
 * A check of saved configuration values by the Worker's Python (`validate_web_config_overrides_json`).
 * `source` says which values: an editor draft's overrides, or a stored request's `config` and
 * `configOverrides`. The port resolves to the Worker's reply, whose values are checked here.
 * @typedef {{ mode: DiagramMode, config: unknown, configOverrides: unknown, managedPaths: readonly string[],
 *   source: 'draft' | 'storedRequest' }} ConfigCheckRequest
 * @typedef {{ result?: { error?: unknown, overrides?: unknown } } | null | undefined} ConfigCheckReply
 * @typedef {(request: ConfigCheckRequest) => ConfigCheckReply | Promise<ConfigCheckReply>} ConfigCheckPort
 */
/** @type {ConfigCheckPort | null} */
let unmanagedConfigOverrideValidator = null;

/** @param {typeof unmanagedConfigOverrideValidator} validator */
export const setUnmanagedConfigOverrideValidator = (validator) => {
  unmanagedConfigOverrideValidator = typeof validator === 'function'
    ? validator
    : null;
};

/**
 * @typedef {{ mode: DiagramMode, config?: unknown, configOverrides?: unknown, source: ConfigCheckRequest['source'] }} UnmanagedOverrideCheck
 *   `mode`: a slice's mode, or the mode of a Session's projection (the projection refuses another value; a
 *   settings-only Session's is checked with its draft). `config` and `configOverrides` are saved values,
 *   unvalidated: the validator service checks them.
 */
/**
 * The overrides the validator service accepts, or its refusal.
 * @param {UnmanagedOverrideCheck} options
 * @returns {Promise<{ overrides: Record<string, unknown> } | { refusal: Error }>}
 */
const checkUnmanagedConfigOverrides = async ({
  mode,
  config = null,
  configOverrides = {},
  source
}) => {
  const managedPaths = managedConfigOverridePathsForMode(mode);
  const overrides = isPlainObject(configOverrides)
    ? cloneJsonData(configOverrides)
    : configOverrides;
  const overridePaths = isPlainObject(overrides) ? Object.keys(overrides) : [];
  if (
    config === null
    && overridePaths.length === 0
  ) {
    return { overrides: {} };
  }
  if (
    config === null
    && source !== 'draft'
    && overridePaths.every((path) => managedPaths.includes(path))
  ) {
    return { overrides: {} };
  }
  if (!unmanagedConfigOverrideValidator) {
    throw new Error('Configuration validation service is unavailable.');
  }
  const result = await unmanagedConfigOverrideValidator({
    mode,
    config,
    configOverrides: overrides,
    managedPaths,
    source
  });
  if (result?.result?.error && typeof result.result.error === 'object') {
    return { refusal: Object.assign(new Error('Configuration validation failed.'), result.result.error) };
  }
  if (typeof result?.result?.error === 'string' && result.result.error.trim()) {
    return { refusal: new Error(result.result.error.trim()) };
  }
  const preserved = result?.result?.overrides;
  if (!isPlainObject(preserved)) {
    throw new Error('Configuration validation returned an invalid preserved-settings result.');
  }
  // The cast: `isPlainObject(preserved)` above (this module's predicate is no type guard).
  return { overrides: cloneJsonData(/** @type {Record<string, unknown>} */ (preserved)) };
};

/**
 * Save and a request: the overrides the validator accepts, else its refusal is thrown.
 * @param {UnmanagedOverrideCheck} options
 * @returns {Promise<Record<string, unknown>>} The overrides the check preserved.
 */
const validateUnmanagedConfigOverrides = async (options) => {
  const checked = await checkUnmanagedConfigOverrides(options);
  if ('refusal' in checked) throw checked.refusal;
  return checked.overrides;
};

/**
 * Load: a stored draft's overrides the validator accepts. When it refuses
 * them, each path is checked alone, and a refused one is dropped with a
 * dialog row (OV-764). A value that is no object is the settings' reader's
 * (`readPlainObject`, listed).
 * @param {UnmanagedOverrideCheck} options
 * @param {string} path Where the overrides are in the file.
 * @param {LoadRepair[]} repairs Written.
 * @returns {Promise<unknown>}
 */
const readStoredUnmanagedConfigOverrides = async (options, path, repairs) => {
  if (!isPlainObject(options.configOverrides)) return options.configOverrides;
  const checked = await checkUnmanagedConfigOverrides(options);
  if (!('refusal' in checked)) return checked.overrides;
  /** @type {Record<string, unknown>} */
  const kept = {};
  for (const [key, value] of Object.entries(/** @type {Record<string, unknown>} */ (options.configOverrides))) {
    const one = await checkUnmanagedConfigOverrides({ ...options, configOverrides: { [key]: value } });
    if ('refusal' in one) repairs.push(loadRepair({ code: 'setting-replaced', path: `${path}["${key}"]`, effect: 'dropped' }));
    else kept[key] = value;
  }
  return validateUnmanagedConfigOverrides({ ...options, configOverrides: kept });
};

/** @param {DrawingState} drawing */
export const applyEditorStateData = (
  drawing,
  editorState = {},
  { normalized: alreadyNormalized = false } = {}
) => {
  const normalized = alreadyNormalized
    ? editorState
    : normalizeEditorStateData(editorState);

  applyEditorArtifactData(normalized);
  applyDrawingEditorData(drawing, normalized);
};

// The shown Result's editor artifacts: its alignment Reset receipt, generated
// Legend inventory and colors, stroke defaults, and feature catalog.
/** @param {Record<string, any>} normalized The result of `normalizeEditorStateData`. */
const applyEditorArtifactData = (normalized) => {
  if (state.similarityAlignmentResetReceipt) {
    state.similarityAlignmentResetReceipt.value = normalized.alignmentResetReceipt ?? null;
  }
  state.originalLegendOrder.value = normalized.legend.originalOrder;
  state.originalLegendColors.value = normalized.legend.originalColors;
  state.originalSvgStroke.value = normalized.originalSvgStroke;
  if (state.featureCatalog) {
    state.featureCatalog.value = admittedFeatureCatalog(normalized.featureCatalog);
  }
};

// The Legend and stroke edits of a drawing, already normalized.
/**
 * @param {DrawingState} drawing
 * @param {Record<string, any>} normalized The result of `normalizeEditorStateData`.
 */
const applyDrawingEditorData = (drawing, normalized) => {
  drawing.legendEntries.value = normalized.legend.entries;
  drawing.deletedLegendEntries.value = normalized.legend.deletedEntries;
  drawing.dormantLegendEntries.value = normalized.legend.dormantEntries;
  replacePlainObject(drawing.legendColorOverrides, normalized.legend.colorOverrides);
  replacePlainObject(drawing.legendStrokeOverrides, normalized.legend.strokeOverrides);
  drawing.addedLegendCaptions.value = new Set(normalized.legend.addedCaptions);
  replacePlainObject(drawing.featureStrokeOverrides, normalized.featureStrokes.overrides);
};

// ---- Session 46: one drawing's slice (`modes.<m>`, PD-OI-086) ----

// A drawing's complete Session 46 slice: every registry row, as Save writes
// it. `selectedFeatureRecordIdx` is the Features-list record of the drawing's
// shown Result (0 for a mode that is not shown).
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 * @param {{ selectedFeatureRecordIdx?: number }} [options]
 */
export const buildModeSliceData = (drawing, mode, { selectedFeatureRecordIdx = 0 } = {}) => {
  const { legend, featureStrokes } = buildEditorStateData(drawing);
  return {
    config: cloneJsonData(buildConfigData(drawing)),
    features: {
      featureOverrides: featureOverridesForState(drawing.featureOverrides),
      featureColorOverrides: cloneJsonData(drawing.featureColorOverrides),
      featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(drawing.featureVisibilityManualRules),
      labelOverrideRows: cloneJsonData(drawing.canonicalLabelOverrideRows.value),
      labelTextBulkOverrides: cloneJsonData(drawing.labelTextBulkOverrides)
    },
    editorState: {
      legend: {
        entries: [...legend.entries, ...legend.dormantEntries.map((entry) => ({ ...entry, dormant: true }))],
        deletedEntries: legend.deletedEntries,
        colorOverrides: legend.colorOverrides,
        strokeOverrides: legend.strokeOverrides,
        addedCaptions: legend.addedCaptions
      },
      featureStrokes
    },
    ui: {
      layoutPreferences: cloneJsonData(drawing.layoutPreferences[mode]),
      canvasPadding: { ...drawing.canvasPadding },
      pendingPaletteName: drawing.pendingPaletteName.value,
      pendingPaletteColors: cloneColors(drawing.pendingPaletteColors.value),
      linearTypographyLinked: Boolean(drawing.linearTypographyLinked.value),
      selectedFeatureRecordIdx
    }
  };
};

const CONFIG_SLICE_ROWS = MODE_SCOPED_SETTINGS.rows.filter((/** @type {{ domain: string }} */ row) => (
  row.domain === 'config' || row.domain.startsWith('config.')
));
// A drawing's settings on Load (plan 4.1): the projection of its mode's
// committed request, then each registry row its slice holds. A slice may
// omit any row; the drawing's defaults fill what neither holds.
/**
 * @param {Record<string, any> | null | undefined} projectedConfig
 * @param {Record<string, any> | null | undefined} sliceConfig
 * @returns {Record<string, any>}
 */
const overlayModeSliceConfig = (projectedConfig, sliceConfig) => {
  /** @type {Record<string, any>} */
  const config = isPlainObject(projectedConfig) ? cloneJsonData(projectedConfig) : {};
  const slice = isPlainObject(sliceConfig) ? sliceConfig : {};
  CONFIG_SLICE_ROWS.forEach((/** @type {{ domain: string, path: string }} */ row) => {
    const container = row.domain === 'config' ? slice : slice[row.domain.slice('config.'.length)];
    if (!isPlainObject(container) || !Object.hasOwn(container, row.path)) return;
    if (row.domain === 'config') {
      config[row.path] = cloneJsonData(container[row.path]);
    } else {
      const domain = row.domain.slice('config.'.length);
      config[domain] = { ...(isPlainObject(config[domain]) ? config[domain] : {}), [row.path]: cloneJsonData(container[row.path]) };
    }
  });
  if (Object.hasOwn(slice, 'colors')) delete config.colorsAreOverrides;
  if (isPlainObject(config.adv)) delete config.adv.losatProgram;
  return config;
};

// Installs one drawing's slice over the projection of its mode's committed
// request (`projectedConfig`, or none); the drawing holds its mode's defaults
// before. `ui` keys apply as Session Load reads them: the pending palette
// only while Instant Preview is off, and the Linear typography link only
// while the two Linear sizes are equal.
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 * @param {Record<string, any>} slice Its `editorState` is read already (`normalizeEditorStateData`, with the rule rows re-keyed).
 * @param {{ projectedConfig?: Record<string, any> | null, resolveTrackPlacements?: boolean, applyCanvasPadding?: boolean }} [options]
 *   `applyCanvasPadding: false` leaves the padding to the caller (Load pads the shown Result after it mounts).
 *   Per-feature color edits are read as the rules are (`savedFeatureColorEdits`).
 */
const applyModeSliceData = (drawing, mode, slice, {
  projectedConfig = null, resolveTrackPlacements = true, applyCanvasPadding = true
} = {}) => {
  const config = overlayModeSliceConfig(projectedConfig, slice.config);
  if (Object.keys(config).length) applyConfigData(drawing, config, { resolveTrackPlacements });
  const features = isPlainObject(slice.features) ? slice.features : {};
  applyDrawingFeatureData(drawing, { ...features, featureColorOverrides: savedFeatureColorEdits(features.featureColorOverrides) });
  const editorState = slice.editorState;
  applyDrawingEditorData(drawing, editorState);
  const ui = isPlainObject(slice.ui) ? slice.ui : {};
  if (isPlainObject(ui.layoutPreferences)) {
    replaceLayoutPreferences(drawing.layoutPreferences, {
      ...cloneJsonData(drawing.layoutPreferences), [mode]: cloneJsonData(ui.layoutPreferences)
    });
  }
  if (applyCanvasPadding && ui.canvasPadding !== undefined) Object.assign(drawing.canvasPadding, readCanvasPadding(ui.canvasPadding).value);
  restorePendingPaletteFromSession(drawing, ui);
  reconcileImportedLinearTypographyLink({ adv: drawing.adv, linked: drawing.linearTypographyLinked, ui });
};

const SESSION_FIELDS_ERROR = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

/**
 * A Session 27-39 (`readSessionDocument` read its band: a Session 31-39
 * holds a `resources` object, and preflight its request object).
 * @template {SessionChecked27 | SessionChecked31 | SessionChecked39} T
 * @param {T} data
 * @param {LoadRepair[]} repairs
 * @returns {T}
 */
const normalizeSessionData = (data, repairs) => {
  const version = data.version;
  return {
    ...data,
    editorState: normalizeEditorStateData(data.editorState, { repairs, sessionVersion: version })
  };
};

/** @template T @param {T} configData @param {boolean} legacy @returns {T} */
const migrateLegacyFeatureRenderingConfig = (configData, legacy) => {
  if (!legacy || !isPlainObject(configData)) return configData;
  const adv = isPlainObject(configData.adv) ? configData.adv : null;
  if (!adv) return configData;
  const features = Array.isArray(adv.features) ? adv.features : null;
  if (features && !features.includes('repeat_region')) return configData;
  const featureShapes = isPlainObject(adv.feature_shapes) ? adv.feature_shapes : {};
  if (Object.prototype.hasOwnProperty.call(featureShapes, 'repeat_region')) return configData;
  return {
    ...configData,
    adv: {
      ...adv,
      feature_shapes: { ...featureShapes, repeat_region: 'rectangle' }
    }
  };
};

/**
 * @param {SessionInput} data
 * @param {'losatCache' | 'losatDerivedCache'} field
 * @returns {unknown[] | null} null: no entry list can be read
 */
const sessionArtifactEntries = (data, field) => {
  const container = data[field];
  if (container === undefined || container === null) return [];
  if (!isPlainObject(container)) return null;
  const entries = Object.prototype.hasOwnProperty.call(container, 'entries') ? container.entries : [];
  return Array.isArray(entries) ? entries : null;
};

// A raw key names one cache entry: an entry with a key an earlier entry has,
// or without one where a key is required, is in `faults`.
/** @param {readonly unknown[]} entries @param {Set<number>} faults @param {{ requireKey?: boolean }} [options] */
const collectInvalidLosatCacheKeys = (entries, faults, { requireKey = false } = {}) => {
  const seen = new Set();
  entries.forEach((entry, index) => {
    const stored = isPlainObject(entry) ? /** @type {{ key?: unknown }} */ (entry).key : '';
    const key = typeof stored === 'string' ? stored : '';
    if (!key) {
      if (requireKey) faults.add(index);
    } else if (seen.has(key)) faults.add(index);
    else seen.add(key);
  });
};

/** A protein identity manifest that holds no protein. */
const EMPTY_PROTEIN_IDENTITY_MANIFEST = Object.freeze({ schema: 2, proteinSets: {}, recordAnalyses: {}, recordInstances: {} });

/**
 * What a Session's LOSAT caches hold that the app cannot read (A13): a cache
 * whose entries are no list (`containers`), the raw and derived entries that
 * cannot be read or repeat a key, and a protein identity manifest that cannot
 * be read (`manifest`), checked as the empty manifest, so the protein entries
 * that need it are faults too. One step per raw entry (Load yields to the page).
 * @typedef {{ containers: Set<'losatCache' | 'losatDerivedCache'>, raw: Set<number>, derived: Set<number>, manifest: boolean }} LosatArtifactFaults
 * @param {SessionInput} data
 * @param {number} sourceSessionVersion
 * @returns {Generator<void, LosatArtifactFaults>}
 */
function* sessionLosatArtifactSteps(data, sourceSessionVersion) {
  /** @type {LosatArtifactFaults} */
  const faults = { containers: new Set(), raw: new Set(), derived: new Set(), manifest: false };
  if (sourceSessionVersion < CURRENT_ARTIFACT_SESSION_MIN_VERSION) return faults;
  const rawRead = sessionArtifactEntries(data, 'losatCache');
  const derivedRead = sessionArtifactEntries(data, 'losatDerivedCache');
  if (rawRead === null) faults.containers.add('losatCache');
  if (derivedRead === null) faults.containers.add('losatDerivedCache');
  const rawEntries = rawRead ?? [];
  const derivedEntries = derivedRead ?? [];
  collectInvalidLosatCacheKeys(rawEntries, faults.raw, { requireKey: true });
  collectInvalidLosatCacheKeys(derivedEntries, faults.derived);

  if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
    recordStructuralMetric('currentSessionPreflightProteinManifestValidationCount');
  }
  let manifest = data.proteinIdentityManifest;
  let identityIndex = buildValidatedProteinIdentityIndex(manifest);
  if (!identityIndex) {
    faults.manifest = true;
    manifest = EMPTY_PROTEIN_IDENTITY_MANIFEST;
    identityIndex = buildValidatedProteinIdentityIndex(manifest);
  }
  try {
    for (const [index, entry] of rawEntries.entries()) {
      if (faults.raw.has(index)) continue;
      const classification = classifyRawLosatCacheEntry(entry);
      if (!['protein-current', 'nucleotide-current'].includes(classification)) {
        faults.raw.add(index);
        continue;
      }
      yield;
      if (classification !== 'protein-current') continue;
      if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
        recordStructuralMetric('currentSessionPreflightProteinRawTextValidationCount');
      }
      if (!validateProteinRawEntryReferences(entry, manifest, { identityIndex })) faults.raw.add(index);
    }
    derivedEntries.forEach((entry, index) => {
      if (!faults.derived.has(index) && !validateDerivedProteinReferences(entry, manifest, { identityIndex })) faults.derived.add(index);
    });
  } finally {
    releaseValidatedProteinIdentityIndex(identityIndex);
  }
  return faults;
};

/** @param {LosatArtifactFaults} faults */
const holdsLosatArtifactFault = (faults) => faults.manifest || faults.containers.size > 0 || faults.raw.size > 0 || faults.derived.size > 0;

/** Save: a Session's LOSAT caches hold no fault. @param {SessionInput} data @param {number} sourceSessionVersion */
export const validateSessionLosatArtifacts = (data, sourceSessionVersion) => {
  const steps = sessionLosatArtifactSteps(data, sourceSessionVersion);
  let step = steps.next();
  while (!step.done) step = steps.next();
  if (holdsLosatArtifactFault(step.value)) throw SESSION_FIELDS_ERROR();
};

/**
 * Load: a Session's LOSAT caches as far as they can be read (A13, E-OD05): a
 * cache whose entries are no list, an entry that cannot be read, and a
 * manifest that cannot be read are dropped (the manifest read as the empty
 * one), each with a row the dialog does not list; the next comparison
 * computes them again. The document is copied where it changed.
 * @template {SessionInput} T
 * @param {T} data
 * @param {number} sourceSessionVersion
 * @param {LoadRepair[]} repairs
 * @returns {Promise<T>}
 */
const readSessionLosatArtifacts = async (data, sourceSessionVersion, repairs) => {
  let deadline = performance.now() + 16;
  const steps = sessionLosatArtifactSteps(data, sourceSessionVersion);
  let step = steps.next();
  while (!step.done) {
    if (performance.now() >= deadline) {
      await new Promise((resolve) => setTimeout(resolve, 0));
      deadline = performance.now() + 16;
    }
    step = steps.next();
  }
  const faults = step.value;
  if (!holdsLosatArtifactFault(faults)) return data;
  /** @type {Record<string, unknown>} */
  const next = { ...data };
  /** @param {string} path */
  const dropped = (path) => repairs.push(loadRepair({ code: 'cache-entry-dropped', path, effect: 'dropped' }));
  faults.containers.forEach((field) => {
    dropped(field);
    delete next[field];
  });
  /** @param {'losatCache' | 'losatDerivedCache'} field @param {Set<number>} indices */
  const withoutEntries = (field, indices) => {
    if (!indices.size || faults.containers.has(field)) return;
    indices.forEach((index) => dropped(`${field}.entries[${index}]`));
    // Not a `containers` fault: the cache is an object whose entries are a list.
    const cache = /** @type {{ entries: unknown[] }} */ (next[field]);
    next[field] = { ...cache, entries: cache.entries.filter((_, index) => !indices.has(index)) };
  };
  withoutEntries('losatCache', faults.raw);
  withoutEntries('losatDerivedCache', faults.derived);
  if (faults.manifest) {
    // A manifest the file does not hold is no loss (E-OD06).
    if (!holdsNothing(data.proteinIdentityManifest)) dropped('proteinIdentityManifest');
    next.proteinIdentityManifest = { ...EMPTY_PROTEIN_IDENTITY_MANIFEST };
  }
  return /** @type {T} */ (next);
};

/**
 * The envelopes of retired LOSAT artifacts a Session keeps, as they were read.
 * @typedef {{ entries?: readonly unknown[] }} LegacyArtifactEnvelope
 */
/**
 * @param {{ legacyRawCandidates?: LegacyArtifactEnvelope | null,
 *   legacyDerivedEvidence?: LegacyArtifactEnvelope | null }} options
 */
export const buildSessionLegacyArtifacts = ({
  legacyRawCandidates,
  legacyDerivedEvidence
}) => {
  /** @type {{ proteinRawCandidates?: LegacyArtifactEnvelope, proteinDerivedEvidence?: LegacyArtifactEnvelope }} */
  const legacyArtifacts = {};
  if (legacyRawCandidates?.entries?.length) {
    legacyArtifacts.proteinRawCandidates = legacyRawCandidates;
  }
  if (legacyDerivedEvidence?.entries?.length) {
    legacyArtifacts.proteinDerivedEvidence = legacyDerivedEvidence;
  }
  return Object.keys(legacyArtifacts).length > 0 ? legacyArtifacts : null;
};

// The Web writers of Sessions 27–33 saved every schema-4 Circular slot row with
// `spacing: null`. With Custom Track Slots off that null is lossless and is
// dropped; the reader table lists any other obsolete slot field (OV-758).
/** @template T @param {T} configData A stored draft. @returns {T} */
const withoutLegacyNullCircularSlotSpacing = (configData) => {
  const adv = /** @type {UnvalidatedFields<'adv'>} */ (configData)?.adv;
  if (!isPlainObject(adv)) return configData;
  // The cast below: `configData` holds `adv`, so it is an object too.
  if (adv.circular_track_slots_enabled
    || adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
    || !Array.isArray(adv.circular_track_slots)) return configData;
  return /** @type {T} */ ({
    .../** @type {Record<string, unknown>} */ (configData),
    adv: {
      ...adv,
      circular_track_slots: adv.circular_track_slots.map((/** @type {unknown} */ slot) => {
        if (!isPlainObject(slot) || slot.spacing !== null) return slot;
        const { spacing: _spacing, ...current } = slot;
        return current;
      })
    }
  });
};


/**
 * A Session's stored drafts (the flat `config` of a Session 27-44 and each
 * Session 46 slice's `config`) with the preserved overrides that are now
 * settings moved into them (`adoptPromotedPreservedOverrides`), read once
 * before anything validates or splits the drafts. A draft whose `adv` is no
 * object takes the font into a new `adv`; the value it drops is listed as
 * the reader of `adv` lists it (OV-823, BA-OD110).
 * @template {object} T
 * @param {T} data
 * @param {LoadRepair[]} repairs
 * @returns {T}
 */
const adoptPromotedPreservedSettings = (data, repairs) => {
  const document = /** @type {{ config?: unknown, modes?: unknown }} */ (data);
  /** @param {unknown} config @param {string} path */
  const adopt = (config, path) => {
    const adoptedConfig = adoptPromotedPreservedOverrides(config);
    // The cast: the adoption changes only a plain object.
    if (adoptedConfig !== config) repairs.push(...valueRepairs(readPlainObject, /** @type {{ adv?: unknown }} */ (config).adv, `${path}.adv`));
    return adoptedConfig;
  };
  /** @type {{ config?: unknown, modes?: unknown }} */
  let adopted = document;
  if (isPlainObject(document.config)) {
    const config = adopt(document.config, 'config');
    if (config !== document.config) adopted = { ...adopted, config };
  }
  if (isPlainObject(document.modes)) {
    const modes = { ...document.modes };
    let changed = false;
    for (const mode of SLICE_MODES) {
      const slice = /** @type {{ config?: unknown }} */ (modes[mode]);
      if (!isPlainObject(slice) || !isPlainObject(slice.config)) continue;
      const config = adopt(slice.config, `modes.${mode}.config`);
      if (config === slice.config) continue;
      modes[mode] = { ...slice, config };
      changed = true;
    }
    if (changed) adopted = { ...adopted, modes };
  }
  return /** @type {T} */ (adopted);
};

/**
 * @param {Record<string, any>} data An unvalidated Session 27-39.
 * @param {number} sourceSessionVersion
 * @returns {Record<string, any>} Its flat draft in the current option names;
 *   Load splits it into mode slices.
 */
const migrateSessionDataToCurrent = (data, sourceSessionVersion) => {
  const readsLegacyOptionValues = sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION;
  const migratedOptions = readsLegacyOptionValues
    ? migratePersistedWebOptionValues(data.config)
    : data.config;
  const circularSlotConfig = migrateImportedCircularTrackSlots(migratedOptions);
  return ({
    ...data,
    version: SESSION_VERSION,
    config: migrateLegacyFeatureRenderingConfig(
      migrateImportedLinearTrackSlots(
        sourceSessionVersion <= 33
          ? withoutLegacyNullCircularSlotSpacing(circularSlotConfig)
          : circularSlotConfig,
        sourceSessionVersion
      ),
      sourceSessionVersion <= 33
    )
  });
};

/**
 * A Session 40 or newer with its drafts' retired field names read as their
 * current names, as an older Session's are (A17g): the drafts' readers read
 * the values.
 * @template {Record<string, unknown>} T
 * @param {T} data
 * @returns {T}
 */
const withCurrentDraftFieldNames = (data) => {
  /** @param {Record<string, unknown>} container A Session, or a slice: its `config` with the current names. */
  const named = (container) => (isPlainObject(container.config)
    ? { ...container, config: migratePersistedWebStateFieldNames(container.config) } : container);
  const modes = isPlainObject(data.modes) ? data.modes : null;
  // The cast: the same Session with its drafts' field names changed.
  return /** @type {T} */ ({
    ...named(data),
    ...(modes ? {
      modes: Object.fromEntries(Object.entries(modes).map(([mode, slice]) => [mode,
        isPlainObject(slice) ? named(slice) : slice]))
    } : {})
  });
};

const LEGACY_CONFIG_KEYS = new Set([
  'form',
  'adv',
  'losat',
  'cliOptions',
  'colors',
  'palette',
  'rules',
  'qualifierPriorityRules',
  'filterMode',
  'whitelist',
  'blacklistText',
  'blastSource',
  'losatProgram',
  'circularConservation'
]);

/**
 * @param {unknown} data A parsed file.
 * @returns {data is Record<string, unknown>}
 */
const isLegacyConfigPayload = (data) =>
  isPlainObject(data) &&
  !Object.prototype.hasOwnProperty.call(data, 'format') &&
  Object.keys(data).some((key) => LEGACY_CONFIG_KEYS.has(key));

// A legacy configuration file, read into the configuration Load applies
// before anything changes (OV-710). Malformed record display settings are
// dropped and listed (BA-D15, OV-569); a track slot that does not read
// refuses the file.
/**
 * @param {Record<string, unknown>} data A legacy configuration (`isLegacyConfigPayload`).
 * @param {DrawingState} drawing The drawing it is applied to.
 * @param {'circular' | 'linear'} mode
 * @returns {Repaired<Record<string, unknown>>}
 */
const readLegacyConfigPayload = (data, drawing, mode) => {
  /** @type {LoadRepair[]} */
  const repairs = [];
  const circularSlotSchema = isPlainObject(data.adv) ? data.adv.circular_track_slots_schema_version : undefined;
  const migratedOptions = circularSlotSchema === CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
    ? data
    : migratePersistedWebOptionValues(data);
  // A shape map that is no object is not used (the reader lists it, and the
  // drawing keeps its shapes), so the repeat-region shape is not added to it.
  const storedShapes = isPlainObject(data.adv) ? data.adv.feature_shapes : undefined;
  /** @type {Record<string, unknown>} */
  const migrated = migrateLegacyFeatureRenderingConfig(
    migrateImportedLinearTrackSlots(
      migrateImportedCircularTrackSlots(migratedOptions)
    ),
    isPlainObject(storedShapes) || isUnsetValue(storedShapes)
  );
  // The reader table reads the track slots' leaves (OV-758, OV-788).
  // What `applyConfigData` reads differently (B37, B43, B48), on a copy of
  // the drawing's settings, which it installs the configuration over.
  /** @type {SettingsContainers} */
  const installed = cloneJsonData({
    form: drawing.form, adv: drawing.adv, losat: drawing.losat, circularConservation: drawing.circularConservation
  });
  repairs.push(...configFieldRepairs(installed, migrated, '', { mode, current: true }));
  // The Legend position of the drawing's mode, read as a Session slice's is.
  if (isPlainObject(migrated.form)) {
    repairs.push(...readLegendPositionOfMode(/** @type {Record<string, unknown>} */ (migrated.form).legend, mode, 'form.legend').repairs);
  }
  return { value: mode === 'linear' ? withBoundLinearRecordsOff(migrated, state.linearSeqs) : migrated, repairs };
};

/**
 * A Linear drawing turns OFF only records it has: a card uid it does not have
 * turned nothing off and is dropped, unlisted.
 * @param {Record<string, unknown>} config
 * @param {readonly { uid?: unknown }[]} linearSeqs
 */
const withBoundLinearRecordsOff = (config, linearSeqs) => {
  if (!Array.isArray(config.recordsOff)) return config;
  const bound = new Set(linearSeqs.map(({ uid }) => uid));
  // Another entry is read (and listed) by the drawing's reader.
  return { ...config, recordsOff: config.recordsOff.filter((uid) => bound.has(uid) || !isRecordDrawKey('linear', uid)) };
};

/** @param {DrawingState} drawing @param {Record<string, unknown>} config What `readLegacyConfigPayload` read. */
const applyLegacyConfigPayload = (drawing, config) => {
  state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(
    drawing, /** @type {UnvalidatedFields<'multi_record_canvas'>} */ (config.form)
  );
  applyConfigData(drawing, config);
  restorePaletteStateAfterConfigImport(drawing);
};

/**
 * @param {DrawingState} drawing
 * @param {UnvalidatedFields<'multi_record_canvas'>} incomingForm
 */
const shouldSuppressCircularMultiRecordDefaults = (drawing, incomingForm) => {
  if (state.mode.value !== 'circular') return false;
  if (!incomingForm || typeof incomingForm !== 'object' || Array.isArray(incomingForm)) return false;
  if (!Object.prototype.hasOwnProperty.call(incomingForm, 'multi_record_canvas')) return false;
  return drawing.form.multi_record_canvas === false && incomingForm.multi_record_canvas === true;
};

/**
 * @param {unknown} stored Unvalidated; anything that is not a plain object counts as empty.
 * @param {unknown} canonical As `stored`.
 * @returns {Record<string, unknown>}
 */
const overlayCanonicalObject = (stored, canonical) => {
  /** @type {Record<string, unknown>} */
  const merged = isPlainObject(stored) ? cloneJsonData(stored) : {};
  if (!isPlainObject(canonical)) return merged;
  Object.entries(canonical).forEach(([key, value]) => {
    if (['__proto__', 'constructor', 'prototype'].includes(key)) return;
    merged[key] = isPlainObject(value)
      ? overlayCanonicalObject(merged[key], value)
      : cloneJsonData(value);
  });
  return merged;
};

/**
 * @param {Record<string, unknown>} projectedConfig
 * @param {unknown} storedConfig Unvalidated.
 * @param {{ hasCanonicalProteinPipeline?: boolean }} [options]
 * @returns {Record<string, unknown>}
 */
const restoreStoredNonCanonicalConfig = (
  projectedConfig,
  storedConfig,
  { hasCanonicalProteinPipeline = false } = {}
) => {
  const restored = cloneJsonData(projectedConfig);
  if (!isPlainObject(storedConfig)) return restored;
  // Keep canonical drawing values authoritative; only supplement state that the
  // canonical request does not currently represent.
  [
    'losat',
    'blastSource',
    'losatProgram',
    'cliOptions',
    'paletteInstantPreviewEnabled',
    'modeProfiles',
    'linearRecordLayout',
    'linearComparisonPlan',
    'webEdits'
  ].forEach((key) => {
    if (hasCanonicalProteinPipeline && key === 'losat') {
      if (!isPlainObject(storedConfig.losat)) return;
      const canonicalLosat = isPlainObject(restored.losat)
        ? cloneJsonData(restored.losat)
        : {};
      restored.losat = overlayCanonicalObject(
        storedConfig.losat,
        canonicalLosat
      );
      return;
    }
    if (
      hasCanonicalProteinPipeline &&
      ['blastSource', 'losatProgram'].includes(key)
    ) return;
    if (Object.prototype.hasOwnProperty.call(storedConfig, key)) {
      restored[key] = cloneJsonData(storedConfig[key]);
    }
  });
  const storedAdv = storedConfig.adv;
  if (isPlainObject(storedAdv) && isPlainObject(restored.adv)) {
    for (const key of [
      'rich_feature_popup',
      'feature_width_circular',
      'depth_width_circular',
      'gc_content_width_circular',
      'gc_content_radius_circular',
      'gc_skew_width_circular',
      'gc_skew_radius_circular'
    ]) {
      if (Object.prototype.hasOwnProperty.call(storedAdv, key)) {
        restored.adv[key] = cloneJsonData(storedAdv[key]);
      }
    }
  }
  const storedConservation = storedConfig?.circularConservation;
  if (!isPlainObject(storedConservation)) {
    return restored;
  }
  const conservation = restored?.circularConservation;
  if (!isPlainObject(conservation)) {
    restored.circularConservation = cloneJsonData(storedConservation);
    return restored;
  }

  ['enabled', 'source', 'losat_program', 'subject_gencode'].forEach((key) => {
    if (Object.prototype.hasOwnProperty.call(storedConservation, key)) {
      conservation[key] = cloneJsonData(storedConservation[key]);
    }
  });
  const series = conservation.series;
  const storedSeries = storedConservation.series;
  if (Array.isArray(series) && Array.isArray(storedSeries)) {
    conservation.series = series.map((/** @type {unknown} */ entry, /** @type {number} */ index) => {
      const storedEntry = storedSeries[index];
      if (!isPlainObject(entry) || !isPlainObject(storedEntry)) return entry;
      if (!Object.prototype.hasOwnProperty.call(storedEntry, 'losat_gencode')) return entry;
      return { ...entry, losat_gencode: cloneJsonData(storedEntry.losat_gencode) };
    });
  }
  return restored;
};

// Version-40 config owns the active controls and the next Generate. The
// canonical projection owns the committed artifact and only fills fields the
// current writer deliberately omitted. Load reads the values with the reader
// table before its dialog (`configFieldRepairs`).
/**
 * @param {{ mode: unknown, projectedConfig: unknown, storedConfig: unknown }} input
 *   Unvalidated; `requireActiveConfigContainers` reads the mode and the stored configuration.
 * @returns {Record<string, unknown>}
 */
export const restoreCurrentWriterActiveConfig = ({
  mode,
  projectedConfig,
  storedConfig
}) => {
  if (!isPlainObject(projectedConfig)) {
    throw new Error('Current session is missing its canonical configuration projection.');
  }
  const activeConfig = requireActiveConfigContainers({ mode, storedConfig });
  const restored = cloneJsonData(projectedConfig);
  /** @type {string[]} */
  const restoredDomains = [];
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS.forEach((domain) => {
    if (!Object.prototype.hasOwnProperty.call(activeConfig, domain)) return;
    restored[domain] = cloneJsonData(activeConfig[domain]);
    restoredDomains.push(domain);
  });
  if (Object.prototype.hasOwnProperty.call(activeConfig, 'colors')) {
    delete restored.colorsAreOverrides;
  }
  if (isPlainObject(restored.adv)) delete restored.adv.losatProgram;

  const projectedRulerLabelFontSize = /** @type {UnvalidatedFields<'ruler_label_font_size'>} */ (projectedConfig.adv)
    ?.ruler_label_font_size;
  if (
    projectedRulerLabelFontSize !== undefined
    && isPlainObject(restored.adv)
    && !Object.prototype.hasOwnProperty.call(restored.adv, 'ruler_label_font_size')
  ) {
    restored.adv.ruler_label_font_size = cloneJsonData(projectedRulerLabelFontSize);
  }

  // Current sessions written before the anchor control became editable omitted
  // this value. Preserve their committed projection only when active config has
  // no value of its own; new sessions store the explicit editor value above.
  const projectedAnchorMode = /** @type {UnvalidatedFields<'blastp', UnvalidatedFields<'collinearAnchorMode'>>} */ (
    projectedConfig.losat)?.blastp?.collinearAnchorMode;
  if (
    projectedAnchorMode !== undefined &&
    isPlainObject(restored.losat) &&
    isPlainObject(restored.losat.blastp) &&
    !Object.prototype.hasOwnProperty.call(restored.losat.blastp, 'collinearAnchorMode')
  ) {
    restored.losat.blastp.collinearAnchorMode = cloneJsonData(projectedAnchorMode);
  }

  recordStructuralMetric('currentWriterActiveConfigRestoreCount', 1, {
    domains: restoredDomains
  });
  recordStructuralMetric('activeConfigCanonicalOverwriteCount', 0, {
    domains: restoredDomains
  });
  return restored;
};

/**
 * @param {SessionInput} data
 * @param {{ adopt?: boolean }} [options]
 */
const validateCurrentWriterFeatureCatalog = async (data, { adopt = false } = {}) => {
  const results = normalizeLogicalResults(
    (Array.isArray(data.results) ? data.results : []).map((/** @type {unknown} */ stored, /** @type {number} */ index) => {
      const result = /** @type {UnvalidatedFields<'name' | 'content'>} */ (stored);
      return { name: result?.name || `Result ${index + 1}`, content: result?.content || '' };
    })
  );
  const catalog = /** @type {UnvalidatedFields<'featureCatalog'>} */ (data.editorState)?.featureCatalog ?? null;
  if (catalog === null) return null;
  return validateFeatureCatalogForImport(catalog, results, {
    adopt,
    mode: /** @type {UnvalidatedFields<'mode'>} */ (data.renderRequest)?.mode || ''
  });
};

/**
 * Load: a committed Result set's feature catalog, checked with its Results
 * by the check Generate and admission refuse with (A9, OV-711): one that does
 * not describe its Results is `dropped`.
 * @param {SessionInput} set
 * @returns {Promise<{ catalog: Awaited<ReturnType<typeof validateCurrentWriterFeatureCatalog>>, dropped: boolean }>}
 */
const readCurrentWriterFeatureCatalog = async (set) => {
  try {
    return { catalog: await validateCurrentWriterFeatureCatalog(set, { adopt: true }), dropped: false };
  } catch (error) {
    if (normalizeCaughtError(error).code !== 'FEATURE_METADATA') throw error;
    return { catalog: null, dropped: true };
  }
};

/**
 * Load: the feature catalog of each committed Result set (A9, OV-711): one
 * that does not describe its Results is dropped. The top-level Results stay
 * without one (`feature-recovery-failed`); the other mode's set, which loads
 * only with its catalog, is dropped whole (`session-part-dropped`).
 * @template {SessionInput} T
 * @param {T} data
 * @param {LoadRepair[]} repairs
 */
const readSessionFeatureCatalogs = async (data, repairs) => {
  const top = await readCurrentWriterFeatureCatalog(data);
  let next = data;
  if (top.dropped) {
    repairs.push(loadRepair({ code: 'feature-recovery-failed', path: 'editorState.featureCatalog', effect: 'dropped' }));
    const editorState = isPlainObject(data.editorState) ? data.editorState : {};
    next = /** @type {T} */ ({ ...data, editorState: { ...editorState, featureCatalog: null } });
  }
  // E1: the other mode's set is validated with its own mode.
  const other = isPlainObject(next.otherModeResult)
    ? await readCurrentWriterFeatureCatalog(/** @type {SessionInput} */ (next.otherModeResult)) : null;
  if (other?.dropped) {
    repairs.push(loadRepair({ code: 'session-part-dropped', path: 'otherModeResult', effect: 'dropped' }));
    const { otherModeResult: _dropped, ...kept } = next;
    next = /** @type {T} */ (kept);
  }
  return { data: next, catalog: top.catalog, otherCatalog: other ? other.catalog : null };
};

/**
 * Load: the alignment reset receipt of each committed Result set, read
 * against the request Load read (RK-2): one that does not bind its plan and
 * sources reads as none, with a row (A19); History, Align and Save refuse it
 * (`validateSimilarityAlignmentResetReceipt`).
 * @template {SessionInput} T
 * @param {T} data
 * @param {LoadRepair[]} repairs
 * @returns {Promise<T>}
 */
const readSessionResetReceipts = async (data, repairs) => {
  /** @param {SessionInput} set @param {string} at @returns {Promise<SessionInput>} */
  const read = async (set, at) => {
    const editorState = isPlainObject(set.editorState) ? set.editorState : {};
    const receipt = editorState.alignmentResetReceipt;
    if (receipt === null || receipt === undefined) return set;
    // Load read the requests and resources before (`readSessionRequests`, `readSessionInputs`).
    const canonical = /** @type {AlignmentResetCanonical} */ ({ renderRequest: set.renderRequest, resources: data.resources });
    const { stale } = await readSimilarityAlignmentResetReceipt(receipt, canonical);
    if (!stale) return set;
    repairs.push(loadRepair({ code: 'session-field-dropped', path: `${at}editorState.alignmentResetReceipt`, effect: 'dropped' }));
    return { ...set, editorState: { ...editorState, alignmentResetReceipt: null } };
  };
  const top = await read(data, '');
  const other = isPlainObject(top.otherModeResult)
    ? await read(/** @type {SessionInput} */ (top.otherModeResult), 'otherModeResult.') : top.otherModeResult;
  return /** @type {T} */ (other === top.otherModeResult ? top : { ...top, otherModeResult: other });
};

/**
 * The stored requests of a Session, read before anything adopts them: the
 * top-level request and `otherModeResult.renderRequest` (OV-670). Returns a
 * new data object when a request changed (the parsed file is never changed),
 * the request schemas the file stored (RK-1), and whether the Session loads as
 * settings only because its request cannot be read (K22, C-OD07): only a
 * Session version that has a settings-only form (42, 44, 46, 47) and holds the
 * settings it keeps does; any other is refused (R4, C-OD25).
 * The Comparisons notice after Load is classified on the stored top-level
 * request, as before the reader (receipt 〔今までの知らせ〕, C-OD14); it is null
 * when Load reads no stored request.
 * @template {object} T
 * @param {T} rawData
 * @param {number} version The Session version Load read.
 * @param {LoadRepair[]} repairs
 * @returns {{ data: T, storedSchemas: StoredRequestSchemas, settingsOnly: boolean,
 *   comparisonClassification: ImportedComparisonIntent | null }}
 */
const readSessionRequests = (rawData, version, repairs) => {
  const top = /** @type {Record<string, unknown>} */ (rawData);
  const other = isPlainObject(top.otherModeResult) ? /** @type {Record<string, unknown>} */ (top.otherModeResult) : null;
  const storedSchemas = {
    renderRequest: /** @type {UnvalidatedFields<'schema'>} */ (top.renderRequest)?.schema,
    otherModeResult: /** @type {UnvalidatedFields<'schema'>} */ (other?.renderRequest)?.schema
  };
  // A missing request is one its reader cannot read (K22, C-OD07; OV-785):
  // no writer omits the key. A null one marks a settings-only Session of a
  // version with that form (a full Session's null request was read as a
  // missing one, `readSessionDocument`); at 31-41 it is one the reader cannot
  // read. A Session before 31 holds no request (no writer wrote one), so Load
  // reads none (OV-862).
  const settingsOnlyForm = hasSettingsOnlyForm(/** @type {SessionEnvelope} */ (rawData));
  if (version < CANONICAL_SESSION_MIN_VERSION || (top.renderRequest === null && settingsOnlyForm)) {
    return { data: rawData, storedSchemas, settingsOnly: false, comparisonClassification: null };
  }
  const resourceIds = new Set(isPlainObject(top.resources) ? Object.keys(/** @type {object} */ (top.resources)) : []);
  /** @param {Record<string, unknown>} set */
  const classify = (set) => classifyImportedComparisonIntent(/** @type {Parameters<typeof classifyImportedComparisonIntent>[0]} */ ({
    renderRequest: set.renderRequest, resources: top.resources
  }));
  /**
   * @param {Record<string, unknown>} set @param {string} path
   * @param {ImportedComparisonIntent} classification
   */
  const read = (set, path, classification) => readStoredRenderRequest(set.renderRequest, {
    featureCatalog: /** @type {UnvalidatedFields<'featureCatalog'>} */ (set.editorState)?.featureCatalog ?? null,
    legacyOrthogroupState: top.orthogroupState ?? null,
    resourceIds,
    path,
    holdComparisons: classification.disposition !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  });
  const comparisonClassification = classify(top);
  const topRead = read(top, 'renderRequest', comparisonClassification);
  if ('kind' in topRead) {
    const envelope = /** @type {SessionEnvelope} */ (rawData);
    repairs.push(loadRepair({ code: 'request-dropped', path: topRead.path, effect: 'dropped' }));
    // Without a settings-only form, an unread request is never committed
    // (C-OD25): the row alone restores nothing, and so does a settings-only
    // reading that keeps no settings (a CLI or Python API Session 42, 44 or
    // 47). R4 refuses both (`restoresNothing`, E-OD22).
    if (!settingsOnlyForm) return { data: rawData, storedSchemas, settingsOnly: false, comparisonClassification: null };
    return {
      data: /** @type {T} */ (settingsOnlySessionCandidate(envelope)), storedSchemas, settingsOnly: true,
      comparisonClassification: null
    };
  }
  repairs.push(...topRead.repairs);
  /** @type {Record<string, unknown>} */
  let data = topRead.value === top.renderRequest ? top : { ...top, renderRequest: topRead.value };
  if (other && Object.hasOwn(other, 'renderRequest')) {
    const otherRead = read(other, 'otherModeResult.renderRequest', classify(other));
    if ('kind' in otherRead) {
      const mode = /** @type {UnvalidatedFields<'mode'>} */ (other.renderRequest)?.mode;
      repairs.push(loadRepair({
        code: 'other-result-set-dropped', path: otherRead.path, effect: 'dropped',
        ...(mode === 'circular' || mode === 'linear' ? { label: mode === 'circular' ? 'Circular' : 'Linear' } : {})
      }));
      const { otherModeResult: _dropped, ...rest } = data;
      data = rest;
    } else {
      repairs.push(...otherRead.repairs);
      if (otherRead.value !== other.renderRequest) data = { ...data, otherModeResult: { ...other, renderRequest: otherRead.value } };
    }
  }
  return { data: /** @type {T} */ (data), storedSchemas, settingsOnly: false, comparisonClassification };
};

/**
 * A schema-2 Gallery Session's request, rebuilt by its migration, read like
 * any stored request: a row here is a defect of that migration.
 * @template {{ renderRequest?: unknown, editorState?: unknown, orthogroupState?: unknown, resources?: unknown }} T
 * @param {T} data
 * @param {LoadRepair[]} repairs
 * @returns {T}
 */
const readRebuiltGalleryRequest = (data, repairs) => {
  const read = readStoredRenderRequest(data.renderRequest, {
    featureCatalog: /** @type {UnvalidatedFields<'featureCatalog'>} */ (data.editorState)?.featureCatalog ?? null,
    legacyOrthogroupState: data.orthogroupState ?? null,
    resourceIds: new Set(isPlainObject(data.resources) ? Object.keys(/** @type {object} */ (data.resources)) : []),
    path: 'renderRequest'
  });
  // The migration rebuilt a request the reader cannot read: nothing is committed unread.
  if ('kind' in read) throw SESSION_FIELDS_ERROR();
  repairs.push(...read.repairs);
  return read.value === data.renderRequest ? data : { ...data, renderRequest: read.value };
};

/**
 * K19: a settings-only Session's shown mode, read before its checks read that
 * mode's settings: a mode of no mode is Circular, with a row; a mode the file
 * does not hold (no `ui`, or a `ui` that holds nothing) is Circular, with none
 * (E-OD06; OV-780).
 * @template {Record<string, unknown>} T
 * @param {T} data @param {LoadRepair[]} repairs
 * @returns {T}
 */
const readSettingsOnlyMode = (data, repairs) => {
  if (!isSettingsOnlySessionDocument(/** @type {SessionEnvelope} */ (/** @type {unknown} */ (data)))) return data;
  const ui = isPlainObject(data.ui) ? data.ui : {};
  if (ui !== data.ui && !holdsNothing(data.ui)) repairs.push(loadRepair({ code: 'session-field-dropped', path: 'ui', effect: 'dropped' }));
  const mode = readSessionUiMode(ui.mode);
  if (mode.lost) repairs.push(...sessionUiFieldRepairs({ mode: ui.mode }, 'ui'));
  return ui === data.ui && ui.mode === mode.value ? data : { ...data, ui: { ...ui, mode: mode.value } };
};

/**
 * @param {SessionEnvelope} stored The document as the file holds it (`restoresNothing`).
 * @param {SessionChecked} sessionData The document `readSessionDocument` read.
 * @param {LoadRepair[]} repairs Collects what Load cannot restore.
 * @param {readonly ResultSetPath[]} uncomposedSets The sets whose Results the composition transform could not compose.
 */
const preflightSessionImport = async (stored, sessionData, repairs, uncomposedSets) => {
  const sourceSessionVersion = sessionData.version;
  // A Session before 47 stored the Legend edits under the shown caption; every
  // reader below reads them by row key (PD-OI-061 revision 5). An edit a later
  // edit of its row overwrote is dropped; a Session 46 drew the later edit, so
  // its drop has a code the dialog does not list (Q1).
  /** @type {{ slice: string, unread: LegendEditsUnread }[]} */
  const unreadLegendEdits = [];
  // A Session before 31 holds no request (OV-862), so no stored table of it is converted.
  const framed = isCanonicalBand(sessionData) ? await convertMainSessionComparisonFrames(sessionData) : sessionData;
  const legendRead = sessionLegendEditsByRowKey(framed, sourceSessionVersion, unreadLegendEdits);
  unreadLegendEdits.forEach(({ slice, unread }) => repairs.push(
    ...legendEditRepairs(unread, `${slice === 'draft' ? '' : `modes.${slice}.`}editorState.legend`)
  ));
  legendRead.outcomes.forEach((outcome) => {
    if (outcome.outcome !== 'collision') return;
    const draft = outcome.slice === 'draft';
    repairs.push(loadRepair({
      code: draft ? 'legend-edit-overwritten' : 'legend-edit-overwritten-session46',
      path: `${draft ? '' : `modes.${outcome.slice}.`}editorState.legend.${outcome.store}`,
      effect: outcome.toKey === null ? 'dropped' : 'changed',
      label: outcome.fromKey,
      used: outcome.rowKey
    }));
  });
  const requestRead = readSessionRequests(legendRead.session, sourceSessionVersion, repairs);
  const requestData = readSettingsOnlyMode(requestRead.data, repairs);
  // R4 decided by the request row of a Session 31-41 (K22): nothing after can restore it.
  if (restoresNothing(stored, repairs)) throw SESSION_FIELDS_ERROR();
  const rawData = await readSessionResetReceipts(
    readSessionAuthorityParts(requestData, sourceSessionVersion, requestRead.storedSchemas, repairs, uncomposedSets), repairs
  );
  const currentSession = isAuthorityBand(rawData);
  /** @type {ReturnType<typeof adoptCurrentSessionDocument> | null} */
  let adoptedSession = null;
  /** @type {ReturnType<typeof adoptCurrentSessionResources> | null} */
  let currentResourceTable = null;
  /** @type {Awaited<ReturnType<typeof validateCurrentWriterFeatureCatalog>> | null} */
  let validatedFeatureCatalog = null;
  /** @type {Awaited<ReturnType<typeof validateCurrentWriterFeatureCatalog>> | null} */
  let otherModeCatalog = null;
  /** @type {ReturnType<typeof normalizeEditorStateData> | null} */
  let normalizedEditorState = null;
  /** @type {SessionChecked} */
  let normalizedData;

  if (currentSession) {
    recordSessionLifecycleEvent('session-authority-validation-start');
    adoptedSession = adoptCurrentSessionDocument(rawData, sourceSessionVersion, requestRead.storedSchemas);
    recordSessionLifecycleEvent('session-authority-validation-end');
    recordSessionLifecycleEvent('feature-catalog-validation-start');
    const catalogRead = await readSessionFeatureCatalogs(rawData, repairs);
    const catalogData = catalogRead.data;
    validatedFeatureCatalog = catalogRead.catalog;
    otherModeCatalog = catalogRead.otherCatalog;
    recordSessionLifecycleEvent('feature-catalog-validation-end');
    recordSessionLifecycleEvent('editor-state-normalization-start');
    normalizedEditorState = normalizeEditorStateData(catalogData.editorState, {
      featureCatalog: validatedFeatureCatalog,
      repairs,
      sessionVersion: sourceSessionVersion
    });
    recordSessionLifecycleEvent('editor-state-normalization-end');
    recordSessionLifecycleEvent('resource-table-adoption-start');
    currentResourceTable = adoptCurrentSessionResources(catalogData.resources);
    recordSessionLifecycleEvent('resource-table-adoption-end');
    // A CLI-written Session 40-41 `config` is no Web draft (OV-269), as in
    // migrate_session_flat_draft.
    const { config: _cliConfig, ...withoutCliConfig } = catalogData;
    normalizedData = withCurrentDraftFieldNames(holdsCliWriterConfig(catalogData) ? withoutCliConfig : catalogData);
  } else {
    validateSessionAuthorityInventory(rawData, sourceSessionVersion, requestRead.storedSchemas, { load: true });
    normalizedData = normalizeSessionData(rawData, repairs);
  }

  // A Gallery Session that stored a schema-2 request (31-33) is rebuilt by
  // its promotion from the request its reader read.
  const promotedData = !currentSession && (
    sourceSessionVersion >= CANONICAL_SESSION_MIN_VERSION &&
    Number(requestRead.storedSchemas.renderRequest) === 2
  )
    ? readRebuiltGalleryRequest(promoteReadGallerySession(normalizedData), repairs)
    : normalizedData;
  if (currentSession) recordSessionLifecycleEvent('losat-artifact-validation-start');
  const cacheData = await readSessionLosatArtifacts(promotedData, sourceSessionVersion, repairs);
  if (currentSession) recordSessionLifecycleEvent('losat-artifact-validation-end');
  const data = adoptPromotedPreservedSettings(currentSession
    ? cacheData
    : migrateSessionDataToCurrent(cacheData, sourceSessionVersion), repairs);
  const storedRequest = /** @type {UnvalidatedFields<'mode'>
    & UnvalidatedFields<'diagramOptions', UnvalidatedFields<'config' | 'configOverrides'>>} */ (data.renderRequest);
  // A Gallery migration rebuilt its request, which is classified as rebuilt.
  const comparisonClassification = promotedData === normalizedData && requestRead.comparisonClassification
    ? requestRead.comparisonClassification
    : classifyImportedComparisonIntent(/** @type {Record<string, any>} */ ({
        renderRequest: data.renderRequest,
        resources: data.resources
      }));
  const projectionRenderRequest = (
    storedRequest?.mode === 'linear'
    && comparisonClassification.disposition
      !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  )
    ? { ...storedRequest, comparisons: [] }
    : storedRequest;
  let currentStoredConfig = sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION
    ? withCurrentLinearLabelVisibility(data.config)
    : data.config;
  if (sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION && isPlainObject(currentStoredConfig)
    && Object.prototype.hasOwnProperty.call(currentStoredConfig, 'recordDisplayDrafts')) {
    currentStoredConfig = {
      ...currentStoredConfig,
      recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(currentStoredConfig.recordDisplayDrafts)
    };
  }
  // A Session 41-44 placement row reached both modes; the migration names
  // its mode (Main rows both), and the Session 46 split moves it there (R2).
  if (sourceSessionVersion < MODE_SCOPED_SESSION_VERSION && isPlainObject(currentStoredConfig)
    && Object.prototype.hasOwnProperty.call(currentStoredConfig, 'featurePlacementOverrides')) {
    currentStoredConfig = {
      ...currentStoredConfig,
      featurePlacementOverrides: migrateSessionFeaturePlacements(currentStoredConfig.featurePlacementOverrides)
    };
  }
  // Session 46 keeps each mode's draft in its slice (PD-OI-086): the committed
  // request's projection reads the slice of the request's mode as a Session
  // 40-44 read its stored draft, and a Session without slices is a CLI or
  // Python writer's.
  const modeScopedSession = sourceSessionVersion >= MODE_SCOPED_SESSION_VERSION;
  const committedSliceConfig = modeScopedSession
    ? /** @type {UnvalidatedFields<'circular' | 'linear', UnvalidatedFields<'config'>>} */ (data.modes)
      ?.[storedRequest?.mode === 'linear' ? 'linear' : 'circular']?.config
    : undefined;
  const hasStoredDraft = modeScopedSession
    ? isPlainObject(committedSliceConfig)
    : Object.prototype.hasOwnProperty.call(data, 'config');
  const runtimeStoredConfig = currentSession && hasStoredDraft
    ? migrateImportedLinearTrackSlots(
        migrateImportedCircularTrackSlots(modeScopedSession ? committedSliceConfig : currentStoredConfig),
        sourceSessionVersion
      )
    : data.config;
  if (currentSession) recordSessionLifecycleEvent('canonical-request-projection-start');
  // The settings-only decision of the request read is passed on explicitly (C-OD07).
  const settingsOnly = requestRead.settingsOnly || isSettingsOnlySessionDocument(data);
  // Kept for the Load stage that projects a Session 39's draft again (OV-1181).
  /** @type {Parameters<typeof projectCanonicalSessionRequest>[0] | null} */
  const canonicalProjectionInput = !settingsOnly && sourceSessionVersion >= CANONICAL_SESSION_MIN_VERSION
    ? {
        renderRequest: projectionRenderRequest,
        resources: data.resources,
        webFiles: data.webFiles,
        legacyFiles: data.files,
        storedConfig: runtimeStoredConfig,
        // A Session without a Web draft (written by the CLI or the Python API)
        // derives its draft from the request (D-02); the comparison plan of a
        // Python re-save is stated the same way as the CLI's (OV-296).
        initializeCliInputs: !hasStoredDraft,
        fileBindings: /** @type {UnvalidatedFields<'fileBindings'>} */ (data.cliInvocation)?.fileBindings,
        linearTrackSlotSchemaVersion: sourceSessionVersion <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
          ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
          : LINEAR_TRACK_SLOT_SCHEMA_VERSION,
        repairLegacyTableRows: sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION,
        sessionResourceTable: currentResourceTable,
        // A Session without a Web draft (CLI, Python) holds its tables only
        // in the request, so Load reads them; a draft holds them already.
        deferResourceContent: currentSession && hasStoredDraft,
        adoptCanonicalPayloads: currentSession,
        // OV-1280: from Session 48 on, a Session without a draft holds each
        // table as its card's file.
        markCardFileRows: !hasStoredDraft && sourceSessionVersion >= FILE_IMPORT_ORIGIN_SESSION_VERSION
      }
    : null;
  /** @type {SessionRequestProjection | null} */
  const canonicalProjection = settingsOnly
    // The cast: Load read the shown mode (`readSettingsOnlyMode`, K19) and Session authority the Web
    // inventory; no step after replaces `ui.mode`.
    ? projectSettingsOnlySession(/** @type {SettingsOnlyProjectionInput} */ (data), currentResourceTable)
    : canonicalProjectionInput && projectCanonicalSessionRequest(canonicalProjectionInput);
  repairs.push(...(canonicalProjection?.legacyTableRepairs || []));
  if (currentSession) recordSessionLifecycleEvent('canonical-request-projection-end');
  let restoredConfig = canonicalProjection
    ? sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
      ? cloneJsonData(canonicalProjection.config)
      : {
          ...restoreStoredNonCanonicalConfig(
            canonicalProjection.config,
            runtimeStoredConfig,
            { hasCanonicalProteinPipeline: Boolean(canonicalProjection.pipelineState) }
          ),
          rules: applySpecificRuleProvenance(
            canonicalProjection.config.rules,
            /** @type {UnvalidatedFields<'rules'>} */ (runtimeStoredConfig)?.rules
          )
        }
    : data.config;
  if (sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION && restoredConfig) {
    restoredConfig = withCurrentLinearLabelVisibility(restoredConfig);
  }
  if (sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION && restoredConfig) {
    const sourceStoredConfig = isPlainObject(normalizedData.config) ? normalizedData.config : {};
    const storedUi = /** @type {UnvalidatedFields<'blastSource'>} */ (normalizedData.ui);
    const forceWebDraft = isPlainObject(sourceStoredConfig.linearRecordLayout)
      || !isPlainObject(sourceStoredConfig.cliOptions);
    const comparisonMigrationConfig = cloneJsonData(restoredConfig);
    if (
      !Object.prototype.hasOwnProperty.call(comparisonMigrationConfig, 'blastSource')
      && storedUi?.blastSource
    ) {
      comparisonMigrationConfig.blastSource = storedUi.blastSource;
    }
    const migratedComparisonDraft = migrateLegacyLinearComparisonDraft({
      config: comparisonMigrationConfig,
      filesData: canonicalProjection?.files || data.files || {},
      forceWebDraft
    });
    restoredConfig = migratedComparisonDraft.config;
    if (canonicalProjection) {
      canonicalProjection.files = migratedComparisonDraft.filesData;
    } else {
      data.files = migratedComparisonDraft.filesData;
    }
  }
  if (canonicalProjection && !modeScopedSession && sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
    && hasStoredDraft) {
    recordSessionLifecycleEvent('current-draft-validation-start');
    restoredConfig = restoreCurrentWriterActiveConfig({
      mode: canonicalProjection.mode,
      projectedConfig: canonicalProjection.config,
      storedConfig: runtimeStoredConfig
    });
    recordSessionLifecycleEvent('current-draft-validation-end');
  }
  // Session 46: each slice is a current-writer draft of its own mode; the
  // committed mode's draft is its slice over its request's projection.
  // `readSessionDocument` read each slice's keys and kinds; the reader table reads its settings.
  if (modeScopedSession) {
    recordSessionLifecycleEvent('current-draft-validation-start');
    restoredConfig = canonicalProjection && !settingsOnly
      ? overlayModeSliceConfig(canonicalProjection.config, runtimeStoredConfig)
      : null;
    if (restoredConfig) {
      // As for Sessions 40-44: the stored draft owns the active controls and the
      // projection fills only what it omits.
      const restoredDomains = CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS
        .filter((domain) => isPlainObject(runtimeStoredConfig) && Object.hasOwn(runtimeStoredConfig, domain));
      recordStructuralMetric('currentWriterActiveConfigRestoreCount', 1, { domains: restoredDomains });
      recordStructuralMetric('activeConfigCanonicalOverwriteCount', 0, { domains: restoredDomains });
    }
    recordSessionLifecycleEvent('current-draft-validation-end');
  }
  if (!canonicalProjection && restoredConfig) {
    hydrateMissingMultiRecordPositionsFromCliInvocation(restoredConfig, /** @type {UnvalidatedFields<'args'>} */ (data.cliInvocation));
  }
  // Saved omission means historical OFF, independently of fresh/reset defaults.
  restoredConfig = {
    ...restoredConfig,
    form: { keep_definition_left_aligned: false, ...restoredConfig?.form }
  };
  const hasCurrentStoredUnmanagedOverrides = currentSession && hasStoredDraft
    && isPlainObject(runtimeStoredConfig)
    && Object.prototype.hasOwnProperty.call(
      runtimeStoredConfig,
      'unmanagedConfigOverrides'
    );
  const unmanagedConfigValidation = canonicalProjection
    ? hasCurrentStoredUnmanagedOverrides
      ? {
          mode: canonicalProjection.mode,
          configOverrides: runtimeStoredConfig.unmanagedConfigOverrides,
          source: /** @type {const} */ ('draft')
        }
      : {
          mode: canonicalProjection.mode,
          config: projectionRenderRequest?.diagramOptions?.config ?? null,
          configOverrides:
            projectionRenderRequest?.diagramOptions?.configOverrides ?? {},
          source: /** @type {const} */ ('storedRequest')
        }
    : null;
  if (currentSession) recordSessionLifecycleEvent('artifact-projection-start');
  const projectionResult = canonicalProjection
    ? (() => {
        const artifactState = projectArtifactState(data);
        if (currentSession) artifactState.editorState = normalizedEditorState;
        return {
          documentMetadata: projectDocumentMetadata(data),
          renderState: {
            mode: canonicalProjection.mode,
            inputType: canonicalProjection.inputType,
            config: restoredConfig,
            semanticFeatureState: canonicalProjection.semanticFeatureState
          },
          editorMetadata: projectWebOnlyEditorMetadata(data),
          artifactState,
          restoredFiles: canonicalProjection.files,
          validatedFeatureCatalog
        };
      })()
    : null;
  if (currentSession) recordSessionLifecycleEvent('artifact-projection-end');
  return {
    data,
    sourceSessionVersion,
    canonicalProjection,
    canonicalProjectionInput,
    restoredConfig,
    projectionResult,
    adoptedCanonicalSession: adoptedSession?.canonical || null,
    currentResourceTable,
    otherModeCatalog,
    comparisonClassification,
    unmanagedConfigValidation,
    settingsOnly
  };
};

// E1: the match-sequence sources one saved Result set (the top-level set or
// `otherModeResult`) needs beyond its catalog's: when the catalog does not
// cover the set's request, or a Session older than 40 has none, the sources
// are read again from the restored files in the set's own mode.
/**
 * @param {{ mode: 'circular' | 'linear', catalog: FeatureCatalog | null, renderRequest: Record<string, any> | undefined,
 *   olderSession: boolean, settingsOnly: boolean, cInputType: string, lInputType: string,
 *   files: Record<string, any>, linearSeqs: Record<string, any>[], circularConservation: Record<string, any> }} set
 * @returns {Promise<{ restored: Record<string, any>[], error: unknown }>}
 */
const restoreLoadedSetSequenceSources = async ({
  mode, catalog, renderRequest, olderSession, settingsOnly, cInputType, lInputType, files, linearSeqs,
  circularConservation
}) => {
  const coverage = catalog
    ? analyzeCatalogSequenceSourceCoverage({
        mode,
        catalogFeatureState: catalog,
        renderRequest,
        comparisonSourceAvailability: mode === 'circular'
          ? resolveCircularComparisonSequenceAvailability({ files, circularConservation })
          : undefined
      })
    : null;
  if (settingsOnly || (!olderSession && coverage?.complete)) return { restored: [], error: null };
  recordStructuralMetric('sourceRecoveryCount');
  try {
    return {
      restored: await buildRestoredMatchSequenceSources({
        mode, cInputType, lInputType, files, linearSeqs, circularConservation
      }),
      error: null
    };
  } catch (sequenceError) {
    console.warn('Session match sequence preparation failed.', normalizeUserFacingError(sequenceError));
    return { restored: [], error: sequenceError };
  }
};

// The stroke and color edits the saved Results of `mode` show (its Session 46
// slice, else the flat draft) and the values the Session kept as Python's
// (`setEditorState`: the set's own editor state). Load records Python's values
// from them on a Result of a Session that saved no such records. `rules` are
// the color rules the mode's drawing loads, as saved: Load names the rule rows
// of a Result an older build saved by their captions (D-19), read as the
// drawing reads them (`savedSpecificRules`).
/**
 * @param {Record<string, any>} data
 * @param {'circular' | 'linear'} mode
 * @param {unknown} setEditorState
 * @param {number} sessionVersion The version of the Session `data`.
 * @param {unknown} [storedRules]
 * @returns {SavedResultEdits}
 */
export const savedResultEdits = (data, mode, setEditorState, sessionVersion, storedRules) => {
  // An entry that is no object is no rule (the rules' reader drops it, listed).
  const rules = Array.isArray(storedRules) ? storedRules.filter(isPlainObject) : [];
  const slice = sessionVersion >= MODE_SCOPED_SESSION_VERSION
    ? data.modes?.[mode]
    : { editorState: data.editorState, features: data.features };
  const edits = normalizeEditorStateData(isPlainObject(slice?.editorState) ? slice.editorState : {});
  const kept = normalizeEditorStateData(isPlainObject(setEditorState) ? /** @type {Record<string, any>} */ (setEditorState) : {});
  const legendRekeys = savedLegendRekeys(rules, kept.legend);
  const legend = rekeyLegendRecords({
    entries: [...edits.legend.entries, ...edits.legend.dormantEntries],
    colorOverrides: edits.legend.colorOverrides,
    strokeOverrides: edits.legend.strokeOverrides,
    originalColors: kept.legend.originalColors,
    originalOrder: kept.legend.originalOrder
  }, legendRekeys);
  return {
    featureColorOverrides: savedFeatureColorEdits(cloneJsonObject(slice?.features?.featureColorOverrides)),
    featureStrokeOverrides: edits.featureStrokes.overrides,
    legendEntries: legend.entries,
    legendColorOverrides: legend.colorOverrides,
    legendStrokeOverrides: legend.strokeOverrides,
    originalLegendColors: legend.originalColors,
    originalLegendOrder: legend.originalOrder,
    rules: savedSpecificRules(rules),
    legendRekeys,
    mode,
    sessionVersion
  };
};

// The rule rows whose key the reading of the saved captions `rules` changes
// (F2), allocated over the generated rows a saved Result set kept (`kept`: its
// inventory and Python's colors; none for a mode without a Result).
/**
 * @param {unknown} rules
 * @param {{ originalOrder: string[], originalColors: Record<string, string> }} kept
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
const savedLegendRekeys = (rules, kept) => {
  /** @type {Map<PythonLegendKey, PythonLegendRow>} */
  const pythonRows = new Map();
  kept.originalOrder.forEach((caption) => {
    const key = /** @type {PythonLegendKey} */ (caption);
    const color = kept.originalColors[caption];
    if (color) pythonRows.set(key, { key, color: /** @type {PythonLegendRow['color']} */ (color) });
  });
  return savedRuleRowRekeys({
    rules: Array.isArray(rules) ? rules : [], pythonRows, features: [], originalLegendOrder: kept.originalOrder
  });
};

// E1: the Results of one saved Result set of a current Session (the top-level
// set or `otherModeResult`), admitted with that set's catalog in its mode.
/**
 * A Result that `legacy.applies` to (no composition metadata: a Session 40
 * Result of main 8228ffab or 7aad9e3e, OV-273) takes `legacy.transform`, the
 * transform of a pre-40 Result; every other Result is admitted as written.
 * @param {Record<string, any>[]} logicalResults
 * @param {{ featureCatalog: FeatureCatalog, mode: 'circular' | 'linear', selectedFeatureTypes: readonly string[] | null | undefined,
 *   savedEdits: SavedResultEdits | null, legacy?: LegacyResultNormalization | null }} set
 */
const admitLoadedSetResults = (logicalResults, { featureCatalog, mode, selectedFeatureTypes, savedEdits, legacy = null }) => {
  const catalogAdmission = admitFeatureCatalog(featureCatalog, logicalResults, { adopt: true, mode });
  return admitCurrentSessionResults(
    createCurrentSessionResultSource(logicalResults, catalogAdmission),
    { mutationPlan: createSavedResultPlan(logicalResults, catalogAdmission, savedEdits, legacy), selectedFeatureTypes }
  );
};

// E1: one saved Result set of a current Session (the top-level set or
// `otherModeResult`) as the artifact slot of its mode, built off-line during
// Load: its Results admitted with its catalog (or the Results Load already
// admitted), its request adopted as that mode's committed Session over the
// shared resource table, and its match-sequence sources.
/**
 * @param {{ renderRequest: Record<string, any>, results: Record<string, any>[], editorState?: unknown,
 *   ui?: unknown, runMetadata?: unknown }} set
 * @param {{ featureCatalog: FeatureCatalog, admittedResults?: Record<string, any>[] | null,
 *   savedEdits?: SavedResultEdits | null, legacy?: LegacyResultNormalization | null,
 *   restoredSequenceSources: Record<string, any>[], resources: Record<string, any>,
 *   resourceTable: any, webFiles: any, retainedBytes: number,
 *   repairs: LoadRepair[] | null }} options `repairs` collects what Load cannot restore of the set's editor state.
 */
const buildLoadedArtifactSlot = (set, {
  featureCatalog, admittedResults = null, savedEdits = null, legacy = null, restoredSequenceSources, resources, resourceTable,
  webFiles, retainedBytes, repairs
}) => {
  const mode = set.renderRequest.mode === 'linear' ? 'linear' : 'circular';
  const results = admittedResults || admitLoadedSetResults(normalizeLogicalResults(set.results.map(
    (/** @type {Record<string, any>} */ result, /** @type {number} */ index) => ({
      name: result?.name || `Result ${index + 1}`,
      content: result?.content || ''
    })
  )), { featureCatalog, mode, selectedFeatureTypes: set.renderRequest.diagramOptions?.selectedFeaturesSet, savedEdits, legacy });
  const features = featureStateFromCatalog(featureCatalog, { mode });
  const editorState = normalizeEditorStateData(
    isPlainObject(set.editorState) ? /** @type {Record<string, any>} */ (set.editorState) : {},
    { featureCatalog, ...(repairs ? { repairs, path: 'otherModeResult.editorState' } : {}) }
  );
  // The set's generated inventory and colors under the keys its rows have after Load (F2).
  const keptLegend = rekeyLegendRecords(editorState.legend, savedEdits?.legendRekeys ?? new Map());
  const ui = isPlainObject(set.ui) ? /** @type {Record<string, any>} */ (set.ui) : {};
  const runMetadata = isPlainObject(set.runMetadata) ? /** @type {Record<string, any>} */ (set.runMetadata) : {};
  const generatedLegend = readLegendPositionOfMode(ui.generatedLegendPosition, mode, 'otherModeResult.ui.generatedLegendPosition');
  repairs?.push(...generatedLegend.repairs, ...sessionUiFieldRepairs(ui, 'otherModeResult.ui'));
  return createArtifactSlot({
    mode,
    values: {
      results,
      selectedResultIndex: Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
        ? Math.min(ui.selectedResultIndex, Math.max(0, results.length - 1)) : 0,
      featureCatalog,
      extractedFeatures: features.extractedFeatures || [],
      biologicalFeatures: features.biologicalFeatures || [],
      featureRecordIds: features.featureRecordIds || [],
      orthogroups: features.orthogroups || [],
      featureOrthogroupIndex: features.featureOrthogroupIndex || new Map(),
      collinearGroups: features.collinearGroups || [],
      trackSlotResolvedGeometry: cloneJsonData(runMetadata.trackSlotGeometry ?? null),
      annotationWarnings: cloneJsonData(runMetadata.annotationWarnings || []),
      featureIdentityNotices: cloneJsonData(runMetadata.featureIdentityNotices || []),
      comparisonWarnings: cloneJsonData(runMetadata.comparisonWarnings || []),
      generatedLegendPosition: generatedLegend.value,
      generatedMode: mode,
      generatedMultiRecordCanvas: truthy(ui.generatedMultiRecordCanvas).value,
      generatedCircularPlotTitlePosition: readGeneratedCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition, ui).value,
      appliedPaletteName: String(ui.appliedPaletteName || 'default'),
      appliedPaletteColors: cloneColors(readPaletteColors(ui.appliedPaletteColors).value),
      similarityAlignmentResetReceipt: editorState.alignmentResetReceipt ?? null,
      originalLegendColors: keptLegend.originalColors,
      originalSvgStroke: editorState.originalSvgStroke
    },
    legendInventory: keptLegend.originalOrder,
    matchSequenceOwner: state.matchSequenceRegistry?.buildSourceOwner?.([
      ...(features.sequenceSources || []), ...restoredSequenceSources
    ]) || null,
    runtimeState: { canonical: {
      committedCanonicalSession: adoptRuntimeCanonicalSession({
        renderRequest: set.renderRequest, resources, webFiles: isPlainObject(webFiles) ? webFiles : {}
      }),
      activeSessionResourceTable: resourceTable
    } },
    retainedBytes
  });
};

// A saved layout owner (current or legacy ui fields) wins; without one, a
// canonical Session takes the layout projected from its committed request
// (`restoredLayoutPreferences`).
/**
 * @param {DrawingState} drawing
 * @param {Record<string, unknown> | null} [ui] A saved `ui` object, unvalidated.
 * @param {{ projected?: StoredLayoutPreferences | null }} [options] The layout projected from the committed request.
 */
const restoreLayoutPreferences = (drawing, ui = {}, { projected = null } = {}) => {
  replaceLayoutPreferences(drawing.layoutPreferences, restoredLayoutPreferences(ui, {
    mode: state.mode.value,
    multiRecord: Boolean(drawing.form.multi_record_canvas),
    projected,
    active: { legend: drawing.form.legend, plotTitlePosition: drawing.adv.plot_title_position }
  }));
};

// The draft fields a label fit row may write: the fields of the config paths
// `gbdraw/labels/fit.py` names for the rows (the label settings and the
// feature overlap switch of "separate overlapping features", LP-D52), never a
// feature-slot input (R10).
/** @type {Readonly<Record<string, (adv: Record<string, unknown>, value: number | string | boolean) => void>>} */
const LABEL_FIT_ADV_WRITERS = Object.freeze({
  outer_label_x_offset: (adv, value) => { adv.outer_label_x_offset = value; },
  outer_label_y_offset: (adv, value) => { adv.outer_label_y_offset = value; },
  label_font_size: (adv, value) => { adv.label_font_size = value; },
  circular_label_spacing: (adv, value) => { adv.circular_label_spacing = value; },
  circular_label_placement: (adv, value) => { adv.circular_label_placement = value; },
  resolve_overlaps: (adv, value) => { adv.resolve_overlaps = value; }
});
/** @type {Readonly<Record<string, (form: Record<string, unknown>, value: string) => void>>} */
const LABEL_FIT_FORM_WRITERS = Object.freeze({
  labels_mode: (form, value) => { form.labels_mode = value; }
});

// Settings a label fit row gives by config path (LP-D43) go into the drawing's
// draft through the request mapping, so the next Generate sends exactly them.
// A path the mode does not manage, or one that is not a label setting, is an
// UNKNOWN_CONFIG_PATH error and writes nothing.
/**
 * @param {DrawingState} drawing
 * @param {Readonly<Record<string, number | string | boolean>>} values
 * @param {'circular' | 'linear'} mode
 */
export const applyLabelFitValues = (drawing, values, mode) => {
  const fields = draftFieldsForConfigValues(values, mode, drawing.adv);
  const form = Object.entries(fields.form);
  const adv = Object.entries(fields.adv);
  if (!form.every(([field]) => Object.hasOwn(LABEL_FIT_FORM_WRITERS, field))
    || !adv.every(([field]) => Object.hasOwn(LABEL_FIT_ADV_WRITERS, field))) {
    throw diagnosticError('INPUT_INVALID', { field: 'configOverrides', reason: 'UNKNOWN_CONFIG_PATH' });
  }
  form.forEach(([field, value]) => LABEL_FIT_FORM_WRITERS[field](drawing.form, value));
  adv.forEach(([field, value]) => LABEL_FIT_ADV_WRITERS[field](drawing.adv, value));
};

// Captured or stored settings (History Undo and Redo, the failed Session Load
// rollback, a settings-only Session) pass `resolveTrackPlacements: false`: the
// track stacks are installed as given, so unset slot sides, lane directions,
// and axis indexes stay unset (R11). Load applies a file's settings with the
// readers it listed their losses with before its dialog (`configFieldRepairs`).
/**
 * @param {DrawingState} drawing
 * @param {Record<string, unknown>} data
 * @param {{ resolveTrackPlacements?: boolean, applyAlignment?: boolean, captured?: boolean }} [options] `captured`:
 *   settings the app wrote (`applyCapturedConfigData`).
 */
const applyConfigData = (drawing, data, { resolveTrackPlacements = true, applyAlignment = true, captured = false } = {}) => {
  const mode = drawing === state.drawings.linear ? 'linear' : 'circular';
  // The execution settings of a Session 44 or older draft or of a request
  // projection are app-level (`state.losatExecution`); the drawing keeps the
  // search settings. A group of no object holds none.
  const losat = readPlainObject(data.losat).value;
  const executionFields = LOSAT_EXECUTION_FIELDS.filter((field) => Object.hasOwn(losat, field));
  if (executionFields.length) {
    applyLosatExecutionData(Object.fromEntries(executionFields.map((field) => [field, losat[field]])));
  }
  // The stored settings' merge and readers (Load lists what they cannot restore first).
  const { arrowRatios } = installStoredSettings(drawing, data, mode, { captured });
  if (arrowRatios) Object.assign(drawing.adv, arrowRatios);
  replacePlainObject(drawing.unmanagedConfigOverrides, cloneJsonData(readPlainObject(data.unmanagedConfigOverrides).value));
  drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length, ...cloneJsonData(readRecordDisplayDrafts(data.recordDisplayDrafts).value));
  // Omitted: every record is drawn.
  drawing.recordsOff.splice(0, drawing.recordsOff.length, ...recordsOffReader(mode)(data.recordsOff).value);
  replacePlainObject(drawing.featurePlacementOverrides, cloneJsonData(featurePlacementsReader(mode)(data.featurePlacementOverrides).value));
  drawing.annotationSets.splice(0, drawing.annotationSets.length, ...readAnnotationSets(data.annotationSets).value);
  const linearLayout = readPlainObject(data.linearRecordLayout).value;
  // Omission takes the fresh default in every Session version.
  drawing.linearRecordLayoutEnabled.value = readLinearRecordLayoutEnabled(linearLayout.enabled).value;
  drawing.linearRecordGap.value = readLinearRecordGap(linearLayout.recordGap).value;
  drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...readLinearRecordRows(linearLayout.rows).value);
  // The record translations and alignment plan are the Linear mode's
  // artifact; only the Linear drawing's settings carry them. A History
  // config step carries neither and leaves them alone (OV-383).
  const ownsAlignment = applyAlignment && drawing === state.drawings.linear;
  if (state.linearRecordTranslations && ownsAlignment) {
    state.linearRecordTranslations.value = cloneJsonData(readList(linearLayout.recordTranslations).value);
  }
  if (state.similarityAlignmentPlan && ownsAlignment) {
    state.similarityAlignmentPlan.value = cloneJsonData(readSimilarityAlignment(linearLayout.similarityAlignment).value);
  }
  replaceLinearComparisonPlan(
    drawing.linearComparisonPlan,
    data.linearComparisonPlan || createDefaultLinearComparisonPlan()
  );
  drawing.importedComparisonIntent.action = readImportedComparisonAction(readPlainObject(data.importedComparisonResolution).value.action).value;
  drawing.form.plot_title = String(drawing.form.plot_title || '');
  // `form.legend` and `adv.plot_title_position` are accessors over
  // `layoutPreferences` (state.js) whose setters normalize a merged value.
  // Writing the resolved value back would only pin an unset Circular
  // multi-record preference, so neither is rewritten here.
  drawing.adv.circular_track_slots_schema_version = CIRCULAR_TRACK_SLOT_SCHEMA_VERSION;
  drawing.adv.circular_track_slots_enabled = drawing.adv.circular_track_slots_enabled === true;
  if (resolveTrackPlacements) {
    const stacks = resolvedTrackStacks(drawing.adv, drawing.form);
    drawing.adv.circular_track_slots_axis_index = stacks.circularAxis;
    drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length, ...stacks.circularSlots);
    drawing.adv.linear_track_slots_schema_version = LINEAR_TRACK_SLOT_SCHEMA_VERSION;
    drawing.adv.linear_track_slots_enabled = drawing.adv.linear_track_slots_enabled === true;
    drawing.adv.linear_track_slots_axis_index = stacks.linearAxis;
    drawing.adv.linear_track_slots.splice(0, drawing.adv.linear_track_slots.length, ...stacks.linearSlots);
  } else {
    // The merge installs the caller's slot objects; copy them so a later edit
    // never writes into a History entry or another caller-owned snapshot.
    ['circular_track_slots', 'linear_track_slots'].forEach((key) => {
      drawing.adv[key].splice(0, drawing.adv[key].length, ...cloneJsonData(drawing.adv[key]));
    });
  }
  if (typeof data.paletteInstantPreviewEnabled === 'boolean') {
    state.paletteInstantPreviewEnabled.value = data.paletteInstantPreviewEnabled;
  }
  const importedPalette = readPalette(data.palette).value;
  if (importedPalette) drawing.selectedPalette.value = importedPalette;
  const storedColors = readPlainObject(data.colors).value;
  if (hasColorEntries(storedColors)) {
    if (data.colorsAreOverrides) {
      const paletteColors = paletteColorsFromDefinitions(drawing.selectedPalette.value) || {};
      drawing.currentColors.value = normalizePaletteColors({
        ...paletteColors,
        ...normalizeColorMap(storedColors)
      });
    } else {
      drawing.currentColors.value = normalizePaletteColors(normalizeColorMap(storedColors));
    }
  } else {
    const paletteColors = paletteColorsFromDefinitions(drawing.selectedPalette.value);
    if (paletteColors) drawing.currentColors.value = paletteColors;
  }

  const storedRules = readColorRules(data.rules);
  if (!storedRules.kept) {
    drawing.manualSpecificRules.splice(0, drawing.manualSpecificRules.length, ...storedRules.value);
    drawing.fileLegendCaptions.value = specificRulesFileCaptions(storedRules.value);
  }
  if (Object.prototype.hasOwnProperty.call(data, 'qualifierPriorityRules')) {
    replaceQualifierPriorityRules(drawing, data.qualifierPriorityRules);
  } else if (Object.prototype.hasOwnProperty.call(data, 'priorityRules')) {
    replaceQualifierPriorityRules(drawing, data.priorityRules);
  }
  if (data.filterMode) drawing.filterMode.value = readFilterMode(data.filterMode).value;
  const storedWhitelist = readWhitelistRules(data.whitelist);
  if (!storedWhitelist.kept) {
    drawing.manualWhitelist.splice(0, drawing.manualWhitelist.length, ...storedWhitelist.value);
  }
  if (data.blacklistText !== undefined) drawing.manualBlacklist.value = readBlacklistText(data.blacklistText).value;
  // Absent before Session 48: no file's colors or words to take back (OV-1280).
  replacePlainObject(drawing.fileImports, readFileImports(data.fileImports).value);
  if (data.losatProgram) {
    drawing.losatProgram.value = readLosatProgram(data.losatProgram).value;
  }
  drawing.circularConservation.labels = String(drawing.circularConservation.labels || '');
  const webEdits = readPlainObject(data.webEdits).value;
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupNameOverrides')) {
    replaceStringMap(drawing.orthogroupNameOverrides, webEdits.orthogroupNameOverrides);
  }
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupDescriptionOverrides')) {
    replaceStringMap(drawing.orthogroupDescriptionOverrides, webEdits.orthogroupDescriptionOverrides);
  }
  // Absent in older Sessions: no dormant names (D-21).
  clearObject(drawing.orthogroupDormantOverrides);
  Object.assign(drawing.orthogroupDormantOverrides, normalizeOrthogroupDormantOverrides(webEdits.orthogroupDormantOverrides));
};

/**
 * Applies captured settings (History Undo and Redo, the Load rollback): they
 * hold only values the app wrote, which `requireCurrentConfigFields` checks.
 * @param {DrawingState} drawing
 * @param {Record<string, unknown>} data
 * @param {{ resolveTrackPlacements?: boolean, applyAlignment?: boolean }} [options]
 */
export const applyCapturedConfigData = (drawing, data, options = {}) => {
  requireCurrentConfigFields(data);
  applyConfigData(drawing, data, { ...options, captured: true });
};

/** @param {DrawingState} drawing */
const restorePaletteStateAfterConfigImport = (drawing) => {
  const draftPaletteName = String(drawing.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const hasPreviewResults = Array.isArray(state.results.value) && state.results.value.length > 0;

  if (
    !hasPreviewResults ||
    state.paletteInstantPreviewEnabled.value ||
    draftPaletteName === String(state.appliedPaletteName.value || '')
  ) {
    state.appliedPaletteName.value = draftPaletteName;
    state.appliedPaletteColors.value = draftColors;
    drawing.pendingPaletteName.value = '';
    drawing.pendingPaletteColors.value = {};
    return;
  }

  drawing.pendingPaletteName.value = draftPaletteName;
  drawing.pendingPaletteColors.value = draftColors;
};

// The palette the shown Result was drawn with (an artifact value), or else
// the drawing's palette.
/**
 * @param {DrawingState} drawing
 * @param {SavedPaletteFields} [ui]
 */
const restoreAppliedPaletteFromSession = (drawing, ui = {}) => {
  const draftPaletteName = String(drawing.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const savedAppliedPaletteName = String(ui.appliedPaletteName || draftPaletteName || 'default');
  const savedAppliedPaletteColors = readPaletteColors(ui.appliedPaletteColors).value ?? draftColors;
  state.appliedPaletteName.value = savedAppliedPaletteName;
  state.appliedPaletteColors.value = normalizePaletteColors(cloneColors(savedAppliedPaletteColors));
};

/** @param {SavedPaletteFields | undefined} ui A saved `ui`. */
const pendingPaletteName = (ui) => String(ui?.pendingPaletteName || '').trim();

// A drawing's pending palette applies only while Instant Preview is off.
/**
 * @param {DrawingState} drawing
 * @param {SavedPaletteFields} [ui]
 */
const restorePendingPaletteFromSession = (drawing, ui = {}) => {
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const savedPendingPaletteName = pendingPaletteName(ui);
  const savedPendingPaletteColors = readPaletteColors(ui.pendingPaletteColors).value ?? draftColors;
  if (!state.paletteInstantPreviewEnabled.value && savedPendingPaletteName) {
    drawing.pendingPaletteName.value = savedPendingPaletteName;
    drawing.pendingPaletteColors.value = normalizePaletteColors(cloneColors(savedPendingPaletteColors));
  } else {
    drawing.pendingPaletteName.value = '';
    drawing.pendingPaletteColors.value = {};
  }
};

/**
 * @param {DrawingState} drawing
 * @param {SavedPaletteFields} [ui]
 */
const restorePaletteStateFromSession = (drawing, ui = {}) => {
  restoreAppliedPaletteFromSession(drawing, ui);
  restorePendingPaletteFromSession(drawing, ui);
};

/**
 * The descriptor each file was saved or loaded as, keyed by that file.
 * @type {WeakMap<object, RawSessionResourceDescriptor>}
 */
const serializedFileDescriptors = new WeakMap();
/**
 * @template {RawSessionResourceDescriptor} D
 * @param {object} file A `File`, or another object whose bytes were read.
 * @param {D} descriptor
 * @returns {D}
 */
const cacheSerializedFileDescriptor = (file, descriptor) => {
  setResourcePayloadOwner(descriptor, file);
  serializedFileDescriptors.set(file, descriptor);
  return descriptor;
};

/**
 * @param {unknown} file A file slot's value: a `File` or a Session resource view. A value without bytes
 *   (a list included) is rejected by `readFileBytes`.
 * @returns {Promise<RawSessionResourceDescriptor | null>}
 */
const serializeFile = async (file) => {
  if (!file) return null;
  // A view: each member is read as an unvalidated value.
  const view = /** @type {UnvalidatedFields<'name' | 'type' | 'lastModified'>} */ (file);
  const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(file));
  if (source?.descriptor) {
    const visibleName = String(view?.name || '').trim();
    const descriptor = visibleName && visibleName !== source.descriptor.name
      ? { ...source.descriptor, name: visibleName }
      : source.descriptor;
    setResourcePayloadOwner(descriptor, file);
    return descriptor;
  }
  // A key is an object (functions included): any other value has no entry.
  const cached = typeof file === 'object' || typeof file === 'function'
    ? serializedFileDescriptors.get(/** @type {object} */ (file))
    : undefined;
  if (cached) return cached;
  const bytes = await readFileBytes(file);
  // `readFileBytes` resolves only for an object (or a function) with `arrayBuffer()`.
  const owner = /** @type {object} */ (file);
  const shared = serializedFileDescriptors.get(owner);
  if (shared) return shared;
  recordStructuralMetric('resourceReencodeCount', 1, {
    resourceName: String(view?.name || 'file')
  });
  recordStructuralMetric('base64EncodeCount', 1, {
    resourceName: String(view?.name || 'file')
  });
  recordStructuralMetric('encodedByteCount', bytes.byteLength, {
    resourceName: String(view?.name || 'file')
  });
  const descriptor = {
    name: view?.name || 'file',
    type: view?.type || '',
    size: bytes.byteLength,
    lastModified: view?.lastModified ?? Date.now(),
    encoding: 'base64',
    data: bytesToBase64(bytes)
  };
  return cacheSerializedFileDescriptor(owner, descriptor);
};

/**
 * What `serializeFileValue` gives for a file slot's value: the descriptor Save stores, or the list of them.
 * @typedef {RawSessionResourceDescriptor | null | SerializedFileValue[]} SerializedFileValue
 */
/**
 * @param {DeserializedFile | undefined} value A file slot's value.
 * @returns {Promise<SerializedFileValue>}
 */
const serializeFileValue = async (value) => (
  Array.isArray(value)
    ? Promise.all(value.map((item) => serializeFileValue(item)))
    : serializeFile(value)
);

/**
 * What `deserializeFile` gives for a stored value: a file slot's value (`FileLike`: a `File` or a Session
 * resource view) or `null`, or a list of them for a stored list.
 * @typedef {FileLike | null | DeserializedFile[]} DeserializedFile
 */
// The `File` constructor converts its name and `type` with ToString and `lastModified` with ToNumber
// (Web IDL), so `String` and `Number` below pass the values it used before; a stored value is JSON or a
// descriptor `serializeFile` made, never a symbol or a BigInt, for which the two conversions differ.
/**
 * @param {unknown} value A stored file value: a descriptor, a Session resource view, or a list of them.
 * @returns {DeserializedFile}
 */
const deserializeFile = (value) => {
  if (Array.isArray(value)) {
    return value.map((/** @type {unknown} */ item) => deserializeFile(item));
  }
  if (isSessionResourceFileView(value)) return value;
  // A view: each member is tested or converted where it is read.
  const entry = /** @type {UnvalidatedFields<'data' | 'name' | 'type' | 'lastModified'>} */ (value);
  if (
    !entry ||
    !Object.prototype.hasOwnProperty.call(entry, 'data') ||
    entry.data === null ||
    entry.data === undefined
  ) return null;
  // A value with an own `data` is an object: a descriptor as the Session stores it.
  const descriptor = /** @type {RawSessionResourceDescriptor} */ (entry);
  if (isEncodedDepthFileEntry(/** @type {{ encoding?: unknown, data?: { schema?: unknown } | null }} */ (entry))) {
    // The test above read `data.schema`; `decodeDepthText` checks every field it reads.
    const text = decodeDepthText(/** @type {DepthFilePayload} */ (entry.data));
    recordStructuralMetric('fileConstructionCount');
    return new File([text], String(entry.name || 'depth.tsv'), {
      type: String(entry.type || 'text/tab-separated-values'),
      lastModified: Number(entry.lastModified ?? Date.now())
    });
  }
  if (typeof entry.data !== 'string') return null;
  const bytes = base64ToBytes(entry.data);
  recordStructuralMetric('base64DecodeCount');
  recordStructuralMetric('decodedByteCount', bytes.byteLength);
  recordStructuralMetric('fileConstructionCount');
  const file = new File([bytes], String(entry.name || 'file'), {
    type: String(entry.type || 'application/octet-stream'),
    lastModified: Number(entry.lastModified ?? Date.now())
  });
  cacheSerializedFileDescriptor(file, descriptor);
  return file;
};

/** @type {{ clearActiveRuntime?: () => void } | null} */
let activePreviewRuntime = null;

/** @param {{ clearActiveRuntime?: () => void } | null | undefined} runtime */
export const setPreviewRuntime = (runtime) => {
  activePreviewRuntime = runtime || null;
};

export const serializeResults = () => {
  if (activePreviewRuntime) {
    return normalizeLogicalResults(state.results.value.map((/** @type {SvgResultEntry} */ res, /** @type {number} */ idx) => ({
      name: res.name || `Result ${idx + 1}`,
      content: res.content
    })));
  }

  const currentSvg = (() => {
    if (!state.svgContainer.value) return null;
    const svg = state.svgContainer.value.querySelector('svg');
    if (!svg) return null;
    return serializeCleanSvg(svg);
  })();

  return normalizeLogicalResults(state.results.value.map((/** @type {SvgResultEntry} */ res, /** @type {number} */ idx) => ({
    name: res.name || `Result ${idx + 1}`,
    content: idx === state.selectedResultIndex.value && currentSvg ? currentSvg : res.content
  })));
};

/**
 * One row of `state.losatCacheInfo`, the pairs a Generate keeps cache entries for: a circular conservation
 * comparison carries `key`, `filename` and `display`, a linear comparison also its edge.
 * @typedef {object} LosatCacheInfoEntry
 * @property {unknown} key The raw cache key: a string from Generate; a Session load stores the saved key (truthy, unchecked).
 * @property {string} filename A file name: Generate's, or the saved one a Session load reads (a default when it is
 *   no text or empty).
 * @property {boolean} display
 * @property {string} [edgeKey]
 * @property {number} [ordinal]
 * @property {string} [queryUid]
 * @property {string} [subjectUid]
 * @property {number} [queryIndex]
 * @property {number} [subjectIndex]
 */
/**
 * The fields of a cache-info row that name a linear pair's comparison edge.
 * @typedef {Omit<LosatCacheInfoEntry, 'key' | 'filename' | 'display'>} LosatCacheInfoIdentity
 */
/** @type {readonly ('edgeKey' | 'queryUid' | 'subjectUid')[]} */
const LOSAT_CACHE_INFO_STRING_FIELDS = ['edgeKey', 'queryUid', 'subjectUid'];
/** @type {readonly ('ordinal' | 'queryIndex' | 'subjectIndex')[]} */
const LOSAT_CACHE_INFO_INTEGER_FIELDS = ['ordinal', 'queryIndex', 'subjectIndex'];
/**
 * @param {NonNullable<UnvalidatedFields<keyof LosatCacheInfoEntry>>} entry A cache-info row, or a saved cache
 *   entry that carries the same fields: unvalidated, each field is tested here. Both callers pass an object.
 */
const losatCacheInfoIdentity = (entry) => {
  /** @type {LosatCacheInfoIdentity} */
  const identity = {};
  LOSAT_CACHE_INFO_STRING_FIELDS.forEach((field) => {
    if (typeof entry?.[field] === 'string' && entry[field].trim()) {
      identity[field] = entry[field];
    }
  });
  LOSAT_CACHE_INFO_INTEGER_FIELDS.forEach((field) => {
    // The cast: `Number.isInteger` in the same statement proves the number.
    if (Number.isInteger(entry?.[field])) identity[field] = /** @type {number} */ (entry[field]);
  });
  return identity;
};

/**
 * @param {DrawingState} drawing
 * @param {NonNullable<UnvalidatedFields<keyof LosatCacheInfoEntry | 'queryRecordInstanceKey'
 *   | 'subjectRecordInstanceKey'>>} entry A saved cache entry that passed `isCurrentRawLosatCacheEntry`; the
 *   fields read here are unvalidated.
 */
const restoredLosatCacheInfoIdentity = (drawing, entry) => {
  const identity = losatCacheInfoIdentity(entry);
  if (state.mode.value !== 'linear') return identity;

  // The search's records: the entry's record instances, or the row's own uids.
  const queryInstanceUid = String(entry?.queryRecordInstanceKey || '').trim();
  const subjectInstanceUid = String(entry?.subjectRecordInstanceKey || '').trim();
  if (Boolean(queryInstanceUid) !== Boolean(subjectInstanceUid)) return {};
  const queryUid = queryInstanceUid || identity.queryUid || '';
  const subjectUid = subjectInstanceUid || identity.subjectUid || '';
  if (!queryUid || !subjectUid || queryUid === subjectUid) return {};

  // Indexes of the committed request: positions in the drawn list.
  const indexByUid = new Map(
    drawnLinearSequences(state.linearSeqs, drawing.recordsOff)
      .map((sequence, index) => [String(sequence?.uid || ''), index])
  );
  const queryIndex = indexByUid.get(queryUid);
  const subjectIndex = indexByUid.get(subjectUid);
  if (queryIndex === undefined || subjectIndex === undefined) return {};

  // The row names the search's edge, which Generate takes in either
  // orientation (OV-987); with no such edge, the search's own pair.
  const resolved = losatEdgeOfPair(drawing.linearComparisonResolution.value.edges, queryIndex, subjectIndex);
  const edge = resolved || {
    edgeKey: linearComparisonEdgeKey(queryUid, subjectUid), queryUid, subjectUid, queryIndex, subjectIndex
  };
  if (
    (identity.edgeKey && identity.edgeKey !== edge.edgeKey)
    || (identity.queryUid && identity.queryUid !== edge.queryUid)
    || (identity.subjectUid && identity.subjectUid !== edge.subjectUid)
  ) return {};
  return {
    ...identity,
    edgeKey: edge.edgeKey,
    queryUid: edge.queryUid,
    subjectUid: edge.subjectUid,
    queryIndex: edge.queryIndex,
    subjectIndex: edge.subjectIndex,
    ...(resolved ? { ordinal: resolved.ordinal } : {})
  };
};

const serializeLosatCache = () => {
  const cacheMap = state.losatCache?.value;
  if (!cacheMap || cacheMap.size === 0) {
    return { entries: [], validatedManifest: null, manifestValidated: false };
  }
  const manifest = state.proteinIdentityManifest.value;
  const rawManifest = toRaw(manifest);
  const info = Array.isArray(state.losatCacheInfo.value) ? state.losatCacheInfo.value : [];
  /** @type {Record<string, unknown>[]} */
  const entries = [];
  const seen = new Set();
  const prevalidatedProteinEntries = new WeakSet();

  /**
   * @param {string} key
   * @param {Record<string, unknown>} cached
   * @param {{ filename?: unknown, display?: unknown }} [infoEntry]
   */
  const buildEntry = (key, cached, infoEntry = {}) => {
    const rawCached = toRaw(cached);
    let serialized = rawCached;
    if (!adoptedLosatCacheValues.has(rawCached)) {
      const { text, ...metadata } = rawCached;
      serialized = {
        ...cloneJsonData(metadata),
        text: String(text ?? '')
      };
    }
    const entry = {
      ...serialized,
      key: String(key),
      filename: String(infoEntry.filename || ''),
      display: Boolean(infoEntry.display),
      ...losatCacheInfoIdentity(infoEntry)
    };
    if (
      classifyRawLosatCacheEntry(entry) === 'protein-current'
      && adoptedLosatCacheValues.get(rawCached) === rawManifest
    ) {
      prevalidatedProteinEntries.add(entry);
    }
    return entry;
  };

  // One entry per raw key: rows that share a key (two rings of one sequence)
  // keep the first row's filename, as the CLI raw cache does.
  info.forEach((/** @type {{ key?: string, filename?: string, display?: boolean } | null} */ entry, /** @type {number} */ idx) => {
    if (!entry || !entry.key || seen.has(entry.key)) return;
    const cached = cacheMap.get(entry.key);
    if (!isCurrentRawLosatCacheEntry(cached)) return;
    entries.push(buildEntry(entry.key, cached, {
      ...entry,
      filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
      display: entry.display !== false
    }));
    seen.add(entry.key);
  });

  cacheMap.forEach((/** @type {Record<string, unknown>} */ value, /** @type {string} */ key) => {
    if (seen.has(key)) return;
    if (!isCurrentRawLosatCacheEntry(value)) return;
    entries.push(buildEntry(key, value));
  });

  // A replaced source can leave earlier bindings in the live cache. Persist
  // only protein evidence that the Session's current manifest can resolve.
  /** @type {ReturnType<typeof buildValidatedProteinIdentityIndex>} */
  let identityIndex = null;
  let reusedValidation = false;
  try {
    return {
      entries: entries.filter((entry) => {
        if (classifyRawLosatCacheEntry(entry) !== 'protein-current') return true;
        if (prevalidatedProteinEntries.has(entry)) {
          reusedValidation = true;
          recordStructuralMetric('sessionSaveProteinRawTextValidationReuseCount');
          return true;
        }
        if (!identityIndex) {
          identityIndex = buildValidatedProteinIdentityIndex(manifest);
        }
        if (!identityIndex) {
          throw new Error('Save Session requires a valid protein identity manifest.');
        }
        recordStructuralMetric('sessionSaveProteinRawTextValidationCount');
        return validateProteinRawEntryReferences(entry, manifest, { identityIndex });
      }),
      validatedManifest: identityIndex || reusedValidation ? manifest : null,
      manifestValidated: Boolean(identityIndex || reusedValidation)
    };
  } finally {
    if (identityIndex) releaseValidatedProteinIdentityIndex(identityIndex);
  }
};

// A saved LOSAT pair cache entry: a current entry with a key, an older
// protein candidate (kept in its envelope), or one Load drops; the next
// comparison computes it again (B46).
/** @param {unknown} entry @returns {'current' | 'legacy' | 'dropped'} */
const savedLosatCacheEntryKind = (entry) => {
  if (classifyRawLosatCacheEntry(entry) === 'protein-legacy') return 'legacy';
  // A current entry is an object (`isCurrentRawLosatCacheEntry`).
  return isCurrentRawLosatCacheEntry(entry) && /** @type {{ key?: unknown }} */ (entry).key ? 'current' : 'dropped';
};

/**
 * @param {DrawingState} drawing
 * @param {unknown} entries
 * @param {unknown} [legacyEnvelope]
 * @param {{ adoptCurrent?: boolean, validatedManifest?: unknown }} [options]
 */
const applyLosatCache = (
  drawing,
  entries,
  legacyEnvelope = null,
  { adoptCurrent = false, validatedManifest = null } = {}
) => {
  const map = new Map();
  /** @type {Record<string, unknown>[]} */
  const info = [];
  /** @type {Record<string, unknown>[]} */
  const legacyEntries = [];

  if (Array.isArray(entries)) {
    entries.forEach((entry, idx) => {
      const kind = savedLosatCacheEntryKind(entry);
      if (kind === 'legacy') legacyEntries.push(entry);
      if (kind !== 'current') return;
      const excludedFields = new Set([
        'key',
        'filename',
        'display',
        ...LOSAT_CACHE_INFO_STRING_FIELDS,
        ...LOSAT_CACHE_INFO_INTEGER_FIELDS
      ]);
      const restored = adoptCurrent
        ? Object.fromEntries(
            Object.entries(entry).filter(([field]) => !excludedFields.has(field))
          )
        : cloneJsonData(entry);
      if (!adoptCurrent) {
        excludedFields.forEach((field) => delete restored[field]);
      } else {
        adoptedLosatCacheValues.set(restored, toRaw(validatedManifest));
      }
      map.set(entry.key, restored);
      if (entry.display === false) return;
      info.push({
        key: entry.key,
        // A saved name that is no text is read as none (OV-981).
        filename: typeof entry.filename === 'string' && entry.filename ? entry.filename : `losat_pair_${idx + 1}.tsv`,
        display: true,
        ...restoredLosatCacheInfoIdentity(drawing, entry)
      });
    });
  }

  state.losatCache.value = map;
  state.losatCacheInfo.value = info;
  const restoredEnvelope = normalizeLegacyProteinCandidateEnvelope(legacyEnvelope);
  const importedEnvelope = createLegacyProteinCandidateEnvelope(legacyEntries);
  state.legacyProteinRawCandidates.value = {
    schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA,
    entries: [...restoredEnvelope.entries, ...importedEnvelope.entries]
  };
};

/** @param {Map<string, unknown> | null | undefined} map */
const pruneLosatDerivedCache = (map) => {
  if (!map || typeof map.delete !== 'function') return;
  while (map.size > LOSAT_DERIVED_CACHE_LIMIT) {
    const oldestKey = map.keys().next().value;
    if (oldestKey === undefined) break;
    map.delete(oldestKey);
  }
};

/**
 * @param {unknown} value Saved evidence, unvalidated. The cast below: `Array.isArray` proves an array; each entry
 *   is tested by the filter.
 * @param {unknown[]} [fallbackEntries]
 */
const normalizeLegacyDerivedEvidence = (value, fallbackEntries = []) => ({
  schema: LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
  entries: [
    ...(
      isPlainObject(value) && value.schema === LEGACY_LOSAT_DERIVED_CACHE_SCHEMA && Array.isArray(value.entries)
        ? /** @type {unknown[]} */ (value.entries)
        : []
    ),
    ...fallbackEntries
  ]
    .filter((entry) => isLosatDerivedCacheEntry(entry) && entry.schema === LEGACY_LOSAT_DERIVED_CACHE_SCHEMA)
    .map((entry) => cloneJsonData(entry))
});

/**
 * The saved derived LOSAT cache as Load reads it: the current entries by key
 * (the newest LOSAT_DERIVED_CACHE_LIMIT), the older evidence entries, and how
 * many entries it drops; the next comparison computes them again (B46).
 * @param {unknown} entries
 * @param {{ adoptCurrent?: boolean }} [options]
 */
const readLosatDerivedCache = (entries, { adoptCurrent = false } = {}) => {
  const map = new Map();
  /** @type {unknown[]} */
  const legacyEntries = [];
  let dropped = 0;
  (Array.isArray(entries) ? entries : []).forEach((entry) => {
    if (!isLosatDerivedCacheEntry(entry)) {
      dropped += 1;
      return;
    }
    if (entry.schema === LEGACY_LOSAT_DERIVED_CACHE_SCHEMA) {
      legacyEntries.push(entry);
      return;
    }
    if (map.has(entry.key)) dropped += 1;
    map.set(entry.key, {
      schema: LOSAT_DERIVED_CACHE_SCHEMA,
      kind: 'derived-losatp-payload',
      idEncoding: 'runtime-handle-v1',
      key: entry.key,
      mode: String(entry.mode || ''),
      payload: adoptCurrent ? entry.payload : cloneJsonData(entry.payload)
    });
  });
  dropped += Math.max(0, map.size - LOSAT_DERIVED_CACHE_LIMIT);
  pruneLosatDerivedCache(map);
  return { map, legacyEntries, dropped };
};

/**
 * @param {unknown} entries
 * @param {unknown} [legacyEvidence]
 * @param {{ adoptCurrent?: boolean }} [options]
 */
const applyLosatDerivedCache = (entries, legacyEvidence = null, { adoptCurrent = false } = {}) => {
  const { map, legacyEntries } = readLosatDerivedCache(entries, { adoptCurrent });
  state.losatDerivedCache.value = map;
  state.legacyProteinDerivedEvidence.value = normalizeLegacyDerivedEvidence(legacyEvidence, legacyEntries);
};

/**
 * @param {unknown} manifest Unvalidated: `validateProteinIdentityManifest` admits it or it is dropped.
 * @param {{ adoptCurrent?: boolean, valid?: boolean }} [options] `valid`: what Load's check found.
 */
const applyProteinIdentityManifest = (manifest, { adoptCurrent = false, valid = validateProteinIdentityManifest(manifest) } = {}) => {
  if (!valid) {
    state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
    adoptedProteinIdentityManifest = null;
    return;
  }
  // `valid`: `validateProteinIdentityManifest` admitted it.
  const admitted = /** @type {ProteinIdentityManifest} */ (manifest);
  state.proteinIdentityManifest.value = adoptCurrent
    ? admitted
    : cloneJsonData(admitted);
  adoptedProteinIdentityManifest = adoptCurrent ? admitted : null;
};

/** @param {unknown} groups Saved similarity groups. @returns {string[]} */
const orthogroupIdsOf = (groups) => (Array.isArray(groups) ? groups : [])
  .map((group) => String(group?.id || '').trim())
  .filter(Boolean);

/**
 * The group name or description edits of the groups a state holds; an edit
 * of another group is dropped (B47).
 * @param {unknown} edits
 * @param {ReadonlySet<string>} groupIds
 */
const orthogroupEditsOf = (edits, groupIds) => {
  const all = Object.entries(cloneStringMap(edits));
  const kept = all.filter(([id]) => groupIds.has(id));
  return { value: Object.fromEntries(kept), dropped: all.length - kept.length };
};

/**
 * @param {DrawingState} drawing
 * @param {Record<string, any>} [orthogroupState]
 * @param {{ legacyRecords?: any, catalogFeatureState?: Record<string, any> | null }} [options]
 */
export const applyOrthogroupStateData = (
  drawing, orthogroupState = {}, { legacyRecords = null, catalogFeatureState = null } = {}
) => {
  const storedGroups = Array.isArray(orthogroupState.groups) ? orthogroupState.groups : [];
  const groups = legacyRecords
    ? migrateLegacyOrthogroupMembers(storedGroups, legacyRecords)
    : storedGroups;
  const groupIds = orthogroupIdsOf(groups);
  const groupIdSet = new Set(groupIds);
  const index = catalogFeatureState?.featureOrthogroupIndex || buildOrthogroupFeatureIndex(groups);

  state.orthogroups.value = groups;
  state.featureOrthogroupIndex.value = index;
  // Current catalog admission has already projected these exact groups onto
  // its features. Reuse that owner rather than deriving the same metadata twice.
  if (!catalogFeatureState) {
    state.extractedFeatures.value = enrichFeaturesWithOrthogroups(
      state.extractedFeatures.value,
      index
    );
    if (state.biologicalFeatures) {
      state.biologicalFeatures.value = enrichFeaturesWithOrthogroups(
        state.biologicalFeatures.value,
        index
      );
    }
  }
  const selectedId = String(orthogroupState.selectedOrthogroupId || '').trim();
  state.selectedOrthogroupId.value = selectedId && groupIdSet.has(selectedId) ? selectedId : (groupIds[0] || '');
  state.selectedOrthogroupAlignmentFeature.value = String(orthogroupState.selectedOrthogroupAlignmentFeature || '').trim();

  replaceStringMap(drawing.orthogroupNameOverrides, orthogroupEditsOf(orthogroupState.orthogroupNameOverrides, groupIdSet).value);
  replaceStringMap(drawing.orthogroupDescriptionOverrides,
    orthogroupEditsOf(orthogroupState.orthogroupDescriptionOverrides, groupIdSet).value);
  if (Object.hasOwn(orthogroupState, 'orthogroupDormantOverrides')) {
    clearObject(drawing.orthogroupDormantOverrides);
    Object.assign(drawing.orthogroupDormantOverrides,
      normalizeOrthogroupDormantOverrides(orthogroupState.orthogroupDormantOverrides));
  }
};

/**
 * The one check that the active mode has its biological inputs; Generate and
 * Save from the draft both use it, so both explain the same missing input.
 * @param {string} [mode] The mode to check; the Generate port passes a drawing's mode as text.
 */
export const assertActiveModeInputs = (mode = state.mode.value, sourceState = state) => {
  const files = sourceState.files || {};
  if (mode === 'circular') {
    const gff = sourceState.cInputType?.value === 'gff';
    if (!(gff ? files.c_gff : files.c_gb)) throw diagnosticError('INPUT_REQUIRED');
    if (gff && !files.c_fasta) throw diagnosticError('FASTA_REQUIRED');
    return;
  }
  const gff = sourceState.lInputType?.value === 'gff';
  const sequences = Array.isArray(sourceState.linearSeqs) ? sourceState.linearSeqs : [];
  sequences.forEach((seq, index) => {
    if (!(gff ? seq?.gff : seq?.gb)) throw diagnosticError('INPUT_REQUIRED', { inputOrdinal: index + 1 });
    if (gff && !seq?.fasta) throw diagnosticError('FASTA_REQUIRED', { inputOrdinal: index + 1 });
  });
  // Only a hand-edited Session turns every record OFF (D-06 asks first).
  if (sequences.length > 0
    && drawnLinearSequences(sequences, sourceState.drawings?.linear?.recordsOff).length === 0) {
    throw diagnosticError('RECORD_SELECTION', { reason: 'NONE_DRAWN' });
  }
};

// `sequences`: the drawn record files, in request order. `cards`: every
// Linear card, OFF cards included; an error names its card number among them.
/**
 * @template {{ uid?: unknown, region_start?: unknown, region_end?: unknown }} S
 * @param {S[]} sequences
 * @param {{ mode?: unknown, status?: unknown, issues?: { code: string, context?: DiagnosticContext }[],
 *   records?: { sourceIndex?: unknown }[] } | null | undefined} catalog
 * @param {{ cards?: { uid: string }[] }} [options]
 * @returns {S[]}
 */
export const materializeLinearRecordFiles = (
  sequences,
  catalog,
  { cards = /** @type {{ uid: string }[]} */ ([]) } = {}
) => {
  const sourceSequences = Array.isArray(sequences) ? sequences : [];
  if (catalog == null) return sourceSequences;
  if (catalog?.mode !== 'linear' || catalog?.status !== 'ready') {
    const issue = Array.isArray(catalog?.issues) ? catalog.issues[0] : null;
    throw issue ? diagnosticError(issue.code, issue.context)
      : diagnosticError('RECORD_SELECTION', { reason: 'DISCOVERY_PENDING' });
  }
  const records = Array.isArray(catalog.records) ? catalog.records : [];
  if (records.length === 0) throw diagnosticError('NO_RECORDS');
  const recordCountBySource = new Map();
  records.forEach((record) => {
    const sourceIndex = Number(record?.sourceIndex);
    recordCountBySource.set(sourceIndex, (recordCountBySource.get(sourceIndex) || 0) + 1);
  });
  sourceSequences.forEach((source, sourceIndex) => {
    const count = recordCountBySource.get(sourceIndex) || 0;
    const inputOrdinal = cards.findIndex((card) => card.uid === source.uid) + 1 || sourceIndex + 1;
    if (count === 0) throw diagnosticError('NO_RECORDS', { inputOrdinal });
    if (count <= 1) return;
    const hasRegion = [source.region_start, source.region_end].some(
      (value) => value !== null && value !== undefined && value !== ''
    );
    if (hasRegion) {
      throw diagnosticError('REGION_INVALID', { inputOrdinal, reason: 'SELECT_RECORD_FOR_REGION' });
    }
  });
  return sourceSequences;
};

/**
 * @param {string} mode
 * @param {Record<string, any>} sourceState The project inputs (`state`).
 * @param {DrawingState} drawing The drawing whose settings choose the inputs.
 * @param {Record<string, any> | null} [comparisonPlanOrOptions]
 */
export const serializeActiveRenderFiles = async (
  mode,
  sourceState,
  drawing,
  comparisonPlanOrOptions = null
) => {
  if (!['circular', 'linear'].includes(mode)) {
    throw new Error(`Unsupported render mode: ${String(mode)}.`);
  }
  const sourceFiles = sourceState.files || {};
  // The request holds the drawn records only (record-selection D-02).
  const normalizedLinearSeqs = mode === 'linear'
    ? normalizeLinearSeqList(drawnLinearSequences(sourceState.linearSeqs, drawing.recordsOff))
    : [];
  const depthRequested = draftDrawsDepth(mode === 'linear' ? 'linear' : 'circular', drawing);
  const serializedLinearSeqs = await Promise.all(
    normalizedLinearSeqs.map(async (seq) => ({
      uid: seq.uid,
      gb: await serializeFile(seq.gb),
      gff: await serializeFile(seq.gff),
      fasta: await serializeFile(seq.fasta),
      depth: depthRequested ? await serializeFileValue(seq.depth) : null,
      losat_gencode: seq.losat_gencode ?? 1,
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
      file_definition: seq.file_definition ?? '',
      file_subtitle: seq.file_subtitle ?? '',
      inferred_definition: seq.inferred_definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }))
  );
  const optionBag = comparisonPlanOrOptions
    && typeof comparisonPlanOrOptions === 'object'
    && (
      Object.prototype.hasOwnProperty.call(comparisonPlanOrOptions, 'comparisonPlan') ||
      Object.prototype.hasOwnProperty.call(comparisonPlanOrOptions, 'linearRecordCatalog')
    )
    ? comparisonPlanOrOptions
    : null;
  const suppliedComparisonPlan = optionBag
    ? optionBag.comparisonPlan
    : comparisonPlanOrOptions;
  const linearSeqs = materializeLinearRecordFiles(
    serializedLinearSeqs,
    optionBag?.linearRecordCatalog ?? null,
    { cards: Array.isArray(sourceState.linearSeqs) ? sourceState.linearSeqs : [] }
  );
  const resolvedComparisonPlan = mode === 'linear'
    ? suppliedComparisonPlan || resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan,
        sequences: normalizeLinearSeqList(sourceState.linearSeqs),
        recordsOff: drawing.recordsOff,
        layout: drawing.linearRecordLayoutEnabled?.value
          ? drawing.linearRecordRows
          : [],
        losatProgram: drawing.losatProgram?.value,
        blastpMode: drawing.losat?.blastp?.mode
      })
    : null;
  const linearComparisons = resolvedComparisonPlan ? await Promise.all(
    // The plan is a resolved one: the caller's snapshot of `resolveLinearComparisonPlan`, or the call above.
    resolvedComparisonPlan.edges
      .filter((/** @type {LinearComparisonResolvedEdge} */ edge) => edge.source === 'upload' && edge.fileActive && edge.file)
      .map(async (/** @type {LinearComparisonResolvedEdge} */ edge) => ({
      id: String(edge.id || ''),
      file: await serializeFile(edge.file)
    }))
  ) : [];
  const linearCanonicalComparisons = mode === 'linear'
    && resolvedComparisonPlan?.mode !== 'none'
    && resolvedComparisonPlan?.edges?.length > 0
    ? await Promise.all(
    (Array.isArray(sourceFiles.linearCanonicalComparisons)
      ? sourceFiles.linearCanonicalComparisons
      : []
    ).map(async (/** @type {ResourceBackedComparison} */ comparison) => (
      isResourceBackedCanonicalComparison(comparison)
        ? {
            ...mapResourceBackedCanonicalComparison(comparison),
            file: await serializeFile(comparison.file)
          }
        : cloneJsonData(comparison)
    ))
    )
    : [];
  const conservationEnabled = mode === 'circular'
    && Boolean(drawing.circularConservation?.enabled);
  const conservationSource = String(drawing.circularConservation?.source || 'upload');
  const includeConservationBlasts = conservationEnabled && (
    conservationSource === 'upload'
    || sourceFiles.c_conservation_blasts_source === 'losat-cache'
  );
  const includeConservationFastas =
    conservationEnabled && conservationSource === 'losat';

  return {
    c_gb: mode === 'circular' ? await serializeFile(sourceFiles.c_gb) : null,
    c_gff: mode === 'circular' ? await serializeFile(sourceFiles.c_gff) : null,
    c_fasta: mode === 'circular' ? await serializeFile(sourceFiles.c_fasta) : null,
    c_depth: mode === 'circular' && depthRequested
      ? await serializeFileValue(sourceFiles.c_depth)
      : null,
    c_conservation_blasts: includeConservationBlasts
      ? await serializeFileValue(sourceFiles.c_conservation_blasts)
      : [],
    c_conservation_blasts_source: includeConservationBlasts
      && sourceFiles.c_conservation_blasts_source === 'losat-cache'
      ? 'losat-cache'
      : null,
    c_conservation_fastas: includeConservationFastas
      ? await serializeFileValue(sourceFiles.c_conservation_fastas)
      : [],
    c_conservation_sequence_sources: includeConservationBlasts
      ? await serializeFileValue(sourceFiles.c_conservation_sequence_sources)
      : [],
    d_color: null,
    t_color: null,
    blacklist: null,
    whitelist: null,
    qualifier_priority: null,
    linearSeqs,
    linearComparisons,
    linearCanonicalComparisons
  };
};

/** A generated comparison as the request describes it. */
const GENERATED_COMPARISON_SHAPE = /** @type {const} */ ({ kind: 'ref', name: 'CanonicalGeneratedProteinComparison' });
/**
 * The saved comparisons of a draft. With `repairs` (Load), an entry that is
 * neither resource-backed nor a generated comparison as the request describes
 * it is dropped and listed (OV-750); Generate passes none.
 * @param {unknown} comparisons
 * @param {{ adoptCanonicalPayloads?: boolean, repairs?: LoadRepair[] | null }} [options]
 */
const deserializeCanonicalComparisons = (
  comparisons,
  { adoptCanonicalPayloads = false, repairs = null } = {}
) => {
  if (!Array.isArray(comparisons)) return [];
  const read = repairs ? comparisons.filter((comparison) => isResourceBackedCanonicalComparison(comparison)
    || shapeViolations(comparison, GENERATED_COMPARISON_SHAPE, CANONICAL_REQUEST_DESCRIPTION.types, '').length === 0)
    : comparisons;
  if (repairs && read.length < comparisons.length) {
    repairs.push(loadRepair({
      code: 'setting-entries-dropped', path: 'files.linearCanonicalComparisons', effect: 'dropped',
      count: comparisons.length - read.length
    }));
  }
  return read.map((comparison) => (
    isResourceBackedCanonicalComparison(comparison)
      ? mapResourceBackedCanonicalComparison(comparison, deserializeFile)
      : adoptCanonicalPayloads ? comparison : cloneJsonData(comparison)
  ));
};

/**
 * @param {{ renderRequest?: unknown, resources?: unknown, webFiles?: unknown }} canonical A Generate candidate:
 *   the request it built or projected, its resources and Web file bindings.
 * @param {{ adoptOwnedRequest?: boolean }} [options]
 */
export const adoptCanonicalRenderArtifacts = (
  canonical,
  { adoptOwnedRequest = false } = {}
) => {
  const ownedCanonical = adoptOwnedRequest
    ? adoptRuntimeCanonicalSession(canonical)
    : canonical;
  const preserveAdoptedResources = isAdoptedCanonicalSession(ownedCanonical);
  const sessionResourceTable = preserveAdoptedResources
    ? adoptCurrentSessionResources(ownedCanonical?.resources)
    : null;
  const projection = projectCanonicalSessionRequest({
    ...ownedCanonical,
    sessionResourceTable,
    deferResourceContent: preserveAdoptedResources,
    adoptCanonicalPayloads: preserveAdoptedResources
  });
  const nextCommittedCanonicalSession = preserveAdoptedResources
    ? adoptRuntimeCanonicalSession(ownedCanonical)
    : cloneCanonicalSession(ownedCanonical);
  /** @type {ReturnType<typeof deserializeCanonicalComparisons> | null} */
  let nextLinearComparisons = null;
  // The Circular conservation settings belong to the Circular drawing.
  const drawing = state.drawings.circular;
  /** @type {{
   *   blasts: any[], fastas: any[], sequenceSources: any[],
   *   projectedConservation: Record<string, any> | undefined, series: Record<string, any>[]
   * } | null} */
  let nextCircularState = null;
  if (projection.mode === 'linear') {
    nextLinearComparisons = deserializeCanonicalComparisons(
      projection.files.linearCanonicalComparisons,
      { adoptCanonicalPayloads: preserveAdoptedResources }
    );
  } else if (projection.files.c_conservation_blasts_source === 'losat-cache') {
    const projectedConservation = projection.config.circularConservation;
    const currentSeries = Array.isArray(drawing.circularConservation.series)
      ? drawing.circularConservation.series.map((/** @type {unknown} */ entry) => cloneJsonData(entry))
      : [];
    const nextSeries = Array.isArray(projectedConservation?.series)
      ? projectedConservation.series.map((entry, index) => ({
          ...entry,
          ...(currentSeries[index] || {}),
          sourceIndex: index,
          fileName: entry.fileName,
          label: entry.label,
          color: entry.color
        }))
      : [];
    nextCircularState = {
      blasts: (projection.files.c_conservation_blasts || [])
        .map((/** @type {unknown} */ entry) => deserializeFile(entry))
        .filter(Boolean),
      fastas: (projection.files.c_conservation_fastas || [])
        .map((/** @type {unknown} */ entry) => deserializeFile(entry)),
      sequenceSources: (projection.files.c_conservation_sequence_sources || [])
        .map((/** @type {unknown} */ entry) => deserializeFile(entry)),
      projectedConservation,
      series: nextSeries
    };
  }

  // Everything that can validate, project, or deserialize finishes before the
  // committed request and its comparison-backed draft resources are replaced.
  activeSessionResourceTable = sessionResourceTable;
  setCommittedCanonicalSession(nextCommittedCanonicalSession);
  if (nextLinearComparisons) {
    state.files.linearCanonicalComparisons = nextLinearComparisons;
  }
  if (nextCircularState) {
    state.files.c_conservation_blasts = nextCircularState.blasts;
    state.files.c_conservation_blasts_source = 'losat-cache';
    state.files.c_conservation_fastas = nextCircularState.fastas;
    state.files.c_conservation_sequence_sources = nextCircularState.sequenceSources;
    if (nextCircularState.projectedConservation) {
      drawing.circularConservation.enabled = true;
      drawing.circularConservation.source = 'losat';
      drawing.circularConservation.reference = nextCircularState.projectedConservation.reference;
      drawing.circularConservation.labels = nextCircularState.series
        .map((entry) => entry.label).join(',');
      drawing.circularConservation.ring_width = nextCircularState.projectedConservation.ring_width;
      drawing.circularConservation.ring_gap = nextCircularState.projectedConservation.ring_gap;
      drawing.circularConservation.series.splice(
        0,
        drawing.circularConservation.series.length,
        ...nextCircularState.series
      );
    }
  }
};

/** @returns {CanonicalRenderRequest | null} */
export const getCommittedCanonicalRenderRequest = () => {
  void committedCanonicalRevision.value;
  return committedCanonicalSession?.renderRequest || null;
};

/**
 * Whether the committed request can show the live value of `requestPath` of
 * its mode's drawing (`liveSettingApplies`); without a committed Session every
 * path can. The app setup gives it to `SettingNote` (`SETTING_NOTE_APPLIES`).
 * @param {string} requestPath
 */
export const committedLiveSettingApplies = (requestPath) => {
  const request = getCommittedCanonicalRenderRequest();
  return !request || liveSettingApplies(request, requestPath, state.drawings[request.mode === 'linear' ? 'linear' : 'circular']);
};

// The cast states what every writer of the variable establishes (see the typedef); it claims no leaf.
/** @returns {CommittedCanonicalSession | null} */
export const getCommittedCanonicalSession = () => /** @type {CommittedCanonicalSession | null} */ (committedCanonicalSession);

// The replacement count of the committed Session: an automatic rerender
// commits what it drew only while it is the one it started from (R1 (c)).
export const readCommittedCanonicalRevision = () => committedCanonicalRevision.value;

// R1 (c): an automatic rerender that applied live setting values commits them
// to the committed Session, which then names what the shown Result draws
// (Save writes them as committed, PD-OI-102). Only those request paths
// change; the resources, the Web file record and the draft stay. The Session
// is replaced, not edited, so a History artifact that holds the previous one
// keeps it.
/** @param {readonly LiveSettingWrite[]} writes */
export const commitLiveSettingValues = (writes) => {
  if (!committedCanonicalSession || writes.length === 0) return;
  const renderRequest = cloneJsonData(committedCanonicalSession.renderRequest);
  applyLiveSettingWrites(renderRequest, writes);
  const next = { ...committedCanonicalSession, renderRequest };
  setCommittedCanonicalSession(isAdoptedCanonicalSession(committedCanonicalSession) ? adoptRuntimeCanonicalSession(next) : next);
};

// The recorded track settings of the shown mode's loaded Result; the other mode has none.
export const getLoadedTrackSettings = () => {
  void committedCanonicalRevision.value;
  return loadedTrackSettings?.mode === state.mode.value ? loadedTrackSettings : null;
};

// The Linear Accession and Length visibility the shown Result was drawn with.
export const getCommittedLinearDefinitionVisibility = () => (
  linearDefinitionVisibilityOf(getCommittedCanonicalRenderRequest())
);

// The record count of a committed source, as the Source recipe reads it before
// it names records by #index (Export Feature Edits TSV).
/**
 * @param {string} resourceId
 * @param {'genbank' | 'fasta'} kind
 */
export const readCommittedResourceRecordCount = (resourceId, kind) => (
  readCanonicalResourceRecordCount(committedCanonicalSession?.resources, resourceId, kind)
);

// `restore(null)` installs an artifact without a committed Session (E1: an
// empty mode slot).
export const canonicalRenderArtifactOwner = Object.freeze({
  capture: () => Object.freeze({ committedCanonicalSession, activeSessionResourceTable }),
  /** @param {{ committedCanonicalSession: any, activeSessionResourceTable: any } | null} snapshot */
  restore: (snapshot) => {
    setCommittedCanonicalSession(snapshot?.committedCanonicalSession ?? null);
    activeSessionResourceTable = snapshot?.activeSessionResourceTable ?? null;
  }
});

/**
 * B-LOAD (OV-1145, D-128): a Session up to 47 stored the protein results of a
 * reversed Linear record in the drawn frame. Once Load shows its stored
 * figure, Python converts them; the bytes replace the stored ones under their
 * IDs in the shown artifact, the kept one and the Linear draft's views of
 * them, before History takes its base. Each resources object and resource
 * table is replaced once. A failure throws, and the Load rolls back (PD-OI-046).
 * @param {number} version
 * @param {Readonly<ArtifactSlot> | null} keptSlot The other mode's artifact.
 * @returns {Promise<Readonly<ArtifactSlot> | null>} The kept artifact holding the converted bytes; null when it is unchanged.
 */
const convertLoadedProteinResultFrames = async (version, keptSlot) => {
  /** @typedef {{ committedCanonicalSession: typeof committedCanonicalSession, activeSessionResourceTable: typeof activeSessionResourceTable }} CanonicalHolder */
  /** @type {CanonicalHolder | null} */
  const kept = keptSlot?.runtimeState?.canonical ?? null;
  /** @type {CanonicalHolder} */
  const shown = { committedCanonicalSession, activeSessionResourceTable };
  const converted = await convertStoredProteinResultFrames(version, [
    shown.committedCanonicalSession, kept?.committedCanonicalSession
  ]);
  if (converted.size === 0) return null;
  /** @type {Map<unknown, Record<string, unknown>>} */
  const resourcesOf = new Map();
  /** @type {Map<unknown, NonNullable<typeof activeSessionResourceTable>>} */
  const tableOf = new Map();
  /** @param {CanonicalHolder} holder @returns {Promise<CanonicalHolder>} */
  const replace = async (holder) => {
    const session = holder.committedCanonicalSession;
    if (!session) return holder;
    if (!resourcesOf.has(session.resources)) {
      resourcesOf.set(session.resources, await withConvertedResources(session.resources, converted));
    }
    const resources = resourcesOf.get(session.resources);
    if (resources === session.resources) return holder;
    const table = holder.activeSessionResourceTable;
    if (table && !tableOf.has(table)) tableOf.set(table, adoptCurrentSessionResources(resources));
    const next = { ...session, resources };
    return {
      committedCanonicalSession: isAdoptedCanonicalSession(session) ? adoptRuntimeCanonicalSession(next) : next,
      activeSessionResourceTable: table ? tableOf.get(table) ?? null : null
    };
  };
  const nextShown = await replace(shown);
  const nextKept = kept ? await replace(kept) : null;
  if (nextShown !== shown) canonicalRenderArtifactOwner.restore(nextShown);
  // The Linear draft reads them through views of the Linear artifact's table.
  const linearTable = [nextShown, nextKept]
    .find((holder) => holder?.committedCanonicalSession?.renderRequest?.mode === 'linear')?.activeSessionResourceTable;
  if (linearTable) {
    /** @type {RequestCanonicalComparison[]} */
    const comparisons = state.files.linearCanonicalComparisons;
    state.files.linearCanonicalComparisons = comparisons.map((comparison) => {
      const view = /** @type {SessionResourceFileView | undefined} */ (comparison?.file);
      const resourceId = /** @type {{ resourceId?: string } | null} */ (getSessionResourceSource(view))?.resourceId;
      return view && resourceId !== undefined && converted.has(resourceId)
        ? { ...comparison, file: createSessionResourceFileView(linearTable, resourceId, view) } : comparison;
    });
  }
  return keptSlot && nextKept && nextKept !== kept
    ? createArtifactSlot({ ...keptSlot, runtimeState: { ...keptSlot.runtimeState, canonical: nextKept } }) : null;
};

/**
 * OV-1181: the same stage for a Session 39, whose Load commits a detached
 * copy and gives the draft detached files. The copy takes the converted
 * bytes, and the Linear draft's comparisons are projected again from them
 * with the input Load projected them with. A failure throws, and the Load
 * rolls back (PD-OI-046).
 * @template {{ renderRequest?: unknown, resources?: unknown }} T
 * @param {number} version
 * @param {T} committed
 * @param {Parameters<typeof projectCanonicalSessionRequest>[0] | null} projectionInput
 * @returns {Promise<T>} `committed` itself when nothing is converted
 */
const convertDetachedProteinResultFrames = async (version, committed, projectionInput) => {
  const converted = await convertStoredProteinResultFrames(version, [/** @type {StoredDrawing} */ (committed)]);
  if (converted.size === 0) return committed;
  const resources = await withConvertedResources(/** @type {Record<string, unknown>} */ (committed.resources), converted);
  const mode = /** @type {UnvalidatedFields<'mode'> | undefined} */ (committed.renderRequest)?.mode;
  if (projectionInput && mode === 'linear') {
    state.files.linearCanonicalComparisons = deserializeCanonicalComparisons(
      projectCanonicalSessionRequest({ ...projectionInput, resources }).files.linearCanonicalComparisons,
      // Load listed the repairs of these entries when it first read them.
      { repairs: [] }
    );
  }
  return { ...committed, resources };
};

/**
 * OV-893 (E-OD30): a Circular comparison series whose BLAST table Load could
 * not read (an empty place of `c_conservation_blasts`) is left out with its
 * files, as Remove leaves a series out (`withoutConservationSeries`), and
 * listed; the series after it keep their own files, and with none left the
 * comparison is off, as after Remove. OV-898 (E-OD35): a place past the
 * draft's series (each place of a draft without series) is named by
 * `labels`, by place, so those labels stand in for its series: the place
 * leaves out its label, and the rings after it keep theirs. Such a place is
 * listed at its table's place when it had a label or an entry Load cannot
 * read; an empty place with neither is no row. Writes the file lists of
 * `files`; returns the draft with the series and labels left.
 * @param {Record<string, unknown>} config The Circular draft Load installs.
 * @param {Record<string, unknown>} files The candidate's file slots (`applyFiles`).
 * @param {{ seriesAt: string, filesAt: string, stored: unknown }} where Where the
 *   draft's series and the input files are in the file, and the stored table
 *   list there (`files.c_conservation_blasts` of a Session 27-30; a dropped
 *   binding of a later one is already listed).
 * @param {LoadRepair[]} repairs
 * @returns {Record<string, unknown>}
 */
const withoutUnreadConservationSeries = (config, files, { seriesAt, filesAt, stored }, repairs) => {
  /** @param {unknown} list @returns {(FileLike | null)[]} */
  const placed = (list) => (Array.isArray(list) ? list : []);
  if (!placed(files.c_conservation_blasts).some((file) => !file)) return config;
  const conservation = isPlainObject(config.circularConservation) ? config.circularConservation : {};
  const series = /** @type {RawConservationSeriesEntry[]} */ (Array.isArray(conservation.series) ? conservation.series : []);
  /** @template {FileLike} T @param {(T | null)[]} list */
  const inSeriesOrder = (list) => orderedOptionalConservationFiles(list, { series });
  // Where each table, in series order, is in the stored list.
  const storedPlaces = inSeriesOrder(placed(files.c_conservation_blasts).map((file, place) => ({ name: file?.name, place })))
    .map((entry) => entry?.place ?? null);
  const storedTables = placed(stored);
  let kept = {
    blasts: inSeriesOrder(placed(files.c_conservation_blasts)),
    fastas: inSeriesOrder(placed(files.c_conservation_fastas)),
    sequenceSources: inSeriesOrder(placed(files.c_conservation_sequence_sources)),
    series: [...series, ...parseConservationLabelText(conservation.labels).slice(series.length).map((label) => ({ label }))],
    labels: '',
    enabled: true
  };
  let keptSeries = series.length;
  for (let index = kept.blasts.length - 1; index >= 0; index -= 1) {
    if (kept.blasts[index]) continue;
    const label = kept.series[index]?.label;
    const place = storedPlaces[index];
    if (index < series.length) {
      keptSeries -= 1;
      repairs.push(loadRepair({
        code: 'comparison-series-dropped', path: `${seriesAt}[${index}]`, effect: 'dropped',
        ...(typeof label === 'string' ? { label } : {})
      }));
    } else if (typeof label === 'string' || (place !== null && !holdsNothing(storedTables[place]))) {
      repairs.push(loadRepair({
        code: 'comparison-series-dropped', path: `${filesAt}.c_conservation_blasts[${place ?? index}]`, effect: 'dropped',
        ...(typeof label === 'string' ? { label } : {})
      }));
    }
    kept = withoutConservationSeries(kept, index);
  }
  files.c_conservation_blasts = kept.blasts;
  files.c_conservation_fastas = kept.fastas;
  files.c_conservation_sequence_sources = kept.sequenceSources;
  // Load turns the comparison off only when no series is left.
  return { ...config, circularConservation: {
    ...conservation, series: kept.series.slice(0, keptSeries), labels: kept.labels, ...(kept.enabled ? {} : { enabled: false })
  } };
};

/**
 * A single file slot's stored value: the first entry `deserializeFile` reads.
 * The other entries of a list, and an entry it cannot read, are dropped and
 * listed at `path` (OV-655, TB23, E-OD07); an entry that holds nothing is no row.
 * @param {unknown} stored
 * @param {string} path
 * @param {LoadRepair[] | null} repairs
 * @returns {FileLike | null}
 */
const readFileSlot = (stored, path, repairs) => {
  const entries = Array.isArray(stored) ? stored : [stored];
  const files = entries.map((entry) => (Array.isArray(entry) ? null : deserializeFile(entry)));
  const kept = files.findIndex(Boolean);
  const dropped = entries.filter((entry, index) => index !== kept && !holdsNothing(entry)).length;
  if (dropped && repairs) repairs.push(loadRepair({ code: 'file-entry-dropped', path, effect: 'dropped', count: dropped }));
  return kept < 0 ? null : /** @type {FileLike} */ (files[kept]);
};

/**
 * A file list's stored value, read by place: an entry `deserializeFile`
 * cannot read keeps its place empty and is listed at `path[index]`; an entry
 * that holds nothing is no row.
 * @param {unknown} stored
 * @param {string} path
 * @param {LoadRepair[] | null} repairs
 */
const readFileList = (stored, path, repairs) => (Array.isArray(stored) ? stored : []).map((entry, index) => {
  const file = deserializeFile(entry);
  if (!file && !holdsNothing(entry)) repairs?.push(loadRepair({ code: 'file-entry-dropped', path: `${path}[${index}]`, effect: 'dropped' }));
  return file;
});

/**
 * `targetDrawing` holds the Linear record layout and comparison plan that the
 * Files reconcile: a drawing, or the Load candidate.
 * @param {Record<string, any> | null} filesData
 * @param {{
 *   adoptCanonicalPayloads?: boolean, resolveRecordInputs?: boolean, targetState?: Record<string, any>,
 *   targetDrawing: Pick<DrawingState, 'linearRecordRows' | 'linearComparisonPlan'>,
 *   repairs?: LoadRepair[] | null, at?: string
 * }} options `repairs`: Load lists what it drops from the saved inputs, at
 *   `at`, where the inputs are in the file (`files`, `webFiles.bindings`).
 */
const applyFiles = (filesData, {
  adoptCanonicalPayloads = false, resolveRecordInputs = true, targetState = state, targetDrawing, repairs = null, at = 'files'
}) => {
  targetState.matchSequenceRegistry?.reset?.();
  targetState.circularRecordList.value = [];
  Object.assign(targetState.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null,
    canonicalRecordIdentities: []
  });
  targetState.files.c_gb = null;
  targetState.files.c_gff = null;
  targetState.files.c_fasta = null;
  targetState.files.c_depth = null;
  targetState.files.c_conservation_blasts = [];
  targetState.files.c_conservation_blasts_source = null;
  targetState.files.c_conservation_fastas = [];
  targetState.files.c_conservation_sequence_sources = [];
  targetState.files.linearCanonicalComparisons = [];
  // The setting-table cards are cleared by Reset Settings, which every Session replace runs first.
  targetState.linearReorderNotice.value = '';

  if (!filesData) {
    targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...normalizeLinearSeqList([]));
    targetDrawing.linearRecordRows.splice(0);
    replaceLinearComparisonPlan(
      targetDrawing.linearComparisonPlan,
      reconcileLinearComparisonPlan(targetDrawing.linearComparisonPlan, targetState.linearSeqs)
    );
    return { collapsedLinearSeqs: false };
  }

  /** @param {string} slot */
  const fileSlot = (slot) => readFileSlot(filesData[slot], `${at}.${slot}`, repairs);
  targetState.files.c_gb = fileSlot('c_gb');
  targetState.files.c_gff = fileSlot('c_gff');
  targetState.files.c_fasta = fileSlot('c_fasta');
  // A Depth slot holds a list.
  targetState.files.c_depth = deserializeFile(filesData.c_depth);
  /** @param {string} slot */
  const fileList = (slot) => readFileList(filesData[slot], `${at}.${slot}`, repairs);
  // A table's place stays until Load leaves out its series (`withoutUnreadConservationSeries`).
  targetState.files.c_conservation_blasts = fileList('c_conservation_blasts');
  targetState.files.c_conservation_blasts_source = filesData.c_conservation_blasts_source === 'losat-cache'
    ? 'losat-cache'
    : null;
  targetState.files.c_conservation_fastas = fileList('c_conservation_fastas');
  targetState.files.c_conservation_sequence_sources = fileList('c_conservation_sequence_sources');
  targetState.files.linearCanonicalComparisons = deserializeCanonicalComparisons(
    filesData.linearCanonicalComparisons,
    { adoptCanonicalPayloads, repairs }
  );
  targetState.files.d_color = fileSlot('d_color');
  targetState.files.t_color = fileSlot('t_color');
  targetState.files.blacklist = fileSlot('blacklist');
  targetState.files.whitelist = fileSlot('whitelist');
  targetState.files.qualifier_priority = fileSlot('qualifier_priority');

  const canonicalCircularRecords = Array.isArray(filesData.circularRecords)
    ? filesData.circularRecords
    : [];
  if (canonicalCircularRecords.length > 0) {
    targetState.circularRecordDiscovery.canonicalRecordIdentities = canonicalCircularRecords
      .map((record, index) => {
        const selector = record?.region?.selector || record?.selector;
        const selectorValue = selector?.kind === 'recordId'
          ? String(selector.value || '')
          : selector?.kind === 'recordIndex'
            ? `#${Number(selector.index) + 1}`
            : `#${index + 1}`;
        return {
          selector: selectorValue,
          record_id: selector?.kind === 'recordId' ? selectorValue : '',
          recordKey: String(record?.recordKey || '')
        };
      });
    Object.assign(targetState.circularRecordDiscovery, {
      status: (targetState.cInputType.value === 'gff'
        ? targetState.files.c_gff && targetState.files.c_fasta
        : targetState.files.c_gb) ? 'deferred' : 'idle',
      error: '',
      inputType: targetState.cInputType.value,
      primaryFile: targetState.cInputType.value === 'gff'
        ? targetState.files.c_gff
        : targetState.files.c_gb,
      pairedFile: targetState.cInputType.value === 'gff'
        ? targetState.files.c_fasta
        : null
    });
  }

  if (Array.isArray(filesData.linearSeqs)) {
    const loadedLinearSeqs = filesData.linearSeqs.map((seq, index) => {
      const seqAt = `${at}.linearSeqs[${index}]`;
      // The genetic code, crop and texts, read once (OV-568).
      const read = readLinearSeqFields(seq, seqAt);
      repairs?.push(...read.repairs);
      return {
        uid: seq.uid,
        gb: readFileSlot(seq.gb, `${seqAt}.gb`, repairs),
        gff: readFileSlot(seq.gff, `${seqAt}.gff`, repairs),
        fasta: readFileSlot(seq.fasta, `${seqAt}.fasta`, repairs),
        depth: deserializeFile(seq.depth),
        ...read.value
      };
    });
    const normalized = normalizeLinearSeqList(loadedLinearSeqs);
    const collapsed = resolveRecordInputs ? collapseEmptyLinearSeqList(loadedLinearSeqs) : loadedLinearSeqs;
    const collapsedLinearSeqs = collapsed.length !== normalized.length;
    targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...collapsed);
    if (resolveRecordInputs) {
      // A saved row, unvalidated (OV-568): the lines below test each field they keep.
      const rowByUid = new Map(targetDrawing.linearRecordRows.map((/** @type {UnvalidatedFields<keyof LinearRecordRow>} */ entry) => [String(entry?.uid || ''), entry]));
      targetDrawing.linearRecordRows.splice(
        0,
        targetDrawing.linearRecordRows.length,
        ...targetState.linearSeqs.map((/** @type {LinearSeq} */ seq, /** @type {number} */ index) => {
          const saved = rowByUid.get(seq.uid);
          const row = Number.isInteger(Number(saved?.row)) && Number(saved.row) > 0
            ? Number(saved.row) : index + 1;
          return {
            uid: seq.uid, row,
            ...(saved?.canonicalCardinality === 'exactly_one'
              ? { canonicalCardinality: 'exactly_one' } : {}),
            ...(Number(saved?.canonicalRow) === row
              && Number.isInteger(saved?.canonicalColumn) && saved.canonicalColumn > 0
              ? { canonicalRow: row, canonicalColumn: saved.canonicalColumn }
              : {})
          };
        })
      );
    }
    const comparisonFiles = new Map(
      /** @type {[string, any][]} */ ((Array.isArray(filesData.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [
          String(comparison?.id || ''),
          deserializeFile(comparison?.file)
        ]))
        .filter(([id]) => id)
    );
    const planWithFiles = normalizeLinearComparisonPlan(targetDrawing.linearComparisonPlan);
    planWithFiles.edges.forEach((edge) => {
      const file = comparisonFiles.get(edge.id) || null;
      // A retained file stays inactive; an edge Load leaves with no file (its
      // binding dropped) is as the pair card's Remove leaves it (OV-899).
      if (file) edge.file = file;
      else setComparisonEdgeFile(edge, null);
    });
    replaceLinearComparisonPlan(
      targetDrawing.linearComparisonPlan,
      reconcileLinearComparisonPlan(planWithFiles, targetState.linearSeqs)
    );
    return { collapsedLinearSeqs };
  }

  targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...normalizeLinearSeqList([]));
  replaceLinearComparisonPlan(
    targetDrawing.linearComparisonPlan,
    reconcileLinearComparisonPlan(targetDrawing.linearComparisonPlan, targetState.linearSeqs)
  );
  return { collapsedLinearSeqs: false };
};

// What the Depth fit of a loaded drawing turns off (B45), read before the
// commit from the settings `applyModeSliceData` installs and the candidate's
// Depth files: the fit keeps every series and every row the app keeps without
// a Depth file (BA-OD94), so only a manual row that names no series of the
// drawing is turned off. It is listed when it was on or the file names its
// series; only the drawing's own mode's stack is listed, as the drawing never
// draws the other.
/**
 * @param {'circular' | 'linear'} mode
 * @param {SettingsContainers} installed The drawing's settings as Load installs them (`configFieldRepairs`).
 * @param {unknown} storedAdv The `adv` the file holds of it.
 * @param {{ at: string, resolveTrackPlacements: boolean, circularFiles: DepthFileSlots | DepthFileSlots[],
 *   linearRows: (FileLike | null)[][] }} options `at`: where `adv` is in the file.
 * @returns {LoadRepair[]}
 */
const depthFitRepairs = (mode, installed, storedAdv, { at, resolveTrackPlacements, circularFiles, linearRows }) => {
  if (!isPlainObject(storedAdv)) return [];
  const stored = /** @type {Record<string, unknown>} */ (storedAdv);
  const installedAdv = installed.adv;
  const installedForm = installed.form;
  const stacks = resolveTrackPlacements ? resolvedTrackStacks(installedAdv, installedForm) : null;
  const fit = fitDepthStateToSources({
    mode,
    adv: stacks ? { ...installedAdv, circular_track_slots: stacks.circularSlots, linear_track_slots: stacks.linearSlots } : installedAdv,
    circularFiles,
    linearRows
  });
  const ownStack = mode === 'linear' ? 'linear_track_slots' : 'circular_track_slots';
  const storedRows = Array.isArray(stored[ownStack]) ? /** @type {unknown[]} */ (stored[ownStack]) : [];
  /** @param {string} id Whether the file's row `id` names its series. */
  const namesSeries = (id) => storedRows.some((row) => {
    // The cast: an object (isPlainObject), whose fields are unchecked.
    const saved = /** @type {{ id?: unknown, params?: { track_index?: unknown } }} */ (isPlainObject(row) ? row : {});
    return String(saved.id || '') === id && isPlainObject(saved.params) && !isUnsetValue(saved.params?.track_index);
  });
  return Object.hasOwn(stored, ownStack) ? fit.disabledRows.filter(({ list, id, enabled }) => list === mode && (enabled || namesSeries(id)))
    .map(({ id }) => loadRepair({ code: 'depth-row-disabled', path: `${at}.${ownStack}`, effect: 'changed', label: id })) : [];
};

// A loaded drawing's Depth series fit its own mode's Depth sources (plan 4.2
// `depth`): each drawing runs this once with its mode. Load listed what the
// fit drops before the commit (`depthFitRepairs`).
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 */
const reconcileDepthTrackStateAfterSessionFiles = (drawing, mode) => {
  const fit = fitDepthStateToSources({
    mode,
    adv: drawing.adv,
    circularFiles: state.files.c_depth,
    linearRows: state.linearSeqs.map((/** @type {LinearSeq} */ seq) => depthFileSlotsFromValue(seq.depth))
  });
  if (fit.linearRowWidth) {
    state.linearSeqs.forEach((/** @type {LinearSeq} */ seq) => {
      seq.depth = padDepthFileSlots(seq.depth, fit.linearRowWidth);
    });
  }
  drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, ...fit.depthTracks);
  drawing.adv.circular_track_slots.splice(
    0,
    drawing.adv.circular_track_slots.length,
    ...applyCircularTrackOrderPlacements(
      fit.circularSlots,
      drawing.adv.nt,
      drawing.form.track_type,
      drawing.adv.circular_track_slots_axis_index
    )
  );
  drawing.adv.linear_track_slots.splice(
    0,
    drawing.adv.linear_track_slots.length,
    ...applyLinearTrackOrderPlacements(
      fit.linearSlots,
      drawing.adv.linear_track_slots_axis_index,
      drawing.adv.nt,
      drawing.form.linear_track_layout
    )
  );
};

/** @param {DrawingState} drawing */
const cloneLiveFileState = (drawing) => ({
  files: {
    ...state.files,
    c_conservation_blasts: Array.isArray(state.files.c_conservation_blasts)
      ? [...state.files.c_conservation_blasts]
      : [],
    c_conservation_fastas: Array.isArray(state.files.c_conservation_fastas)
      ? [...state.files.c_conservation_fastas]
      : [],
    c_conservation_sequence_sources: Array.isArray(state.files.c_conservation_sequence_sources)
      ? [...state.files.c_conservation_sequence_sources]
      : [],
    linearCanonicalComparisons: (
      Array.isArray(state.files.linearCanonicalComparisons)
        ? state.files.linearCanonicalComparisons
        : []
    ).map((/** @type {ResourceBackedComparison} */ comparison) => (
      isResourceBackedCanonicalComparison(comparison)
        ? mapResourceBackedCanonicalComparison(comparison)
        : cloneJsonData(comparison)
    ))
  },
  circularRecordList: cloneJsonData(state.circularRecordList.value),
  circularRecordDiscovery: { ...state.circularRecordDiscovery },
  linearSeqs: state.linearSeqs.map((/** @type {LinearSeq} */ seq) => ({
    ...seq,
    depth: Array.isArray(seq.depth) ? [...seq.depth] : seq.depth
  })),
  linearRecordRows: drawing.linearRecordRows.map((/** @type {LinearRecordRow} */ entry) => ({ ...entry })),
  linearComparisonPlan: {
    mode: drawing.linearComparisonPlan.mode,
    defaultSource: drawing.linearComparisonPlan.defaultSource,
    edges: drawing.linearComparisonPlan.edges.map((/** @type {LinearComparisonEdge} */ edge) => ({ ...edge }))
  }
});

/**
 * @param {DrawingState} drawing
 * @param {{
 *   files: Record<string, unknown>, circularRecordList: unknown, circularRecordDiscovery: object,
 *   linearSeqs: LinearSeq[], linearRecordRows: LinearRecordRow[], linearComparisonPlan: LinearComparisonPlan
 * }} snapshot The project inputs `cloneLiveFileState` captured, or those a Load installs.
 */
const restoreLiveFileState = (drawing, snapshot) => {
  state.matchSequenceRegistry?.reset?.();
  Object.keys(state.files).forEach((key) => {
    state.files[key] = snapshot.files[key] ?? null;
  });
  state.circularRecordList.value = cloneJsonData(snapshot.circularRecordList);
  Object.assign(state.circularRecordDiscovery, snapshot.circularRecordDiscovery);
  state.linearSeqs.splice(0, state.linearSeqs.length, ...snapshot.linearSeqs);
  drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...snapshot.linearRecordRows);
  replaceLinearComparisonPlan(drawing.linearComparisonPlan, snapshot.linearComparisonPlan);
};

const captureSessionImportTransientState = () => ({
  semanticFileWatchersSuppressed: Boolean(
    state.semanticFileWatchersSuppressed.value
  ),
  skipCaptureBaseConfig: Boolean(state.skipCaptureBaseConfig.value),
  suppressCircularMultiRecordDefaults: Boolean(
    state.suppressCircularMultiRecordDefaults.value
  ),
  linearReorderNotice: state.linearReorderNotice.value,
  linearFoldedSources: [...state.linearFoldedSources],
  rightDrawer: captureRightDrawerState(state),
  showCanvasControls: Boolean(state.showCanvasControls.value),
  isPanning: Boolean(state.isPanning.value),
  panStart: cloneJsonData(state.panStart),
  selectedAnnotation: state.selectedAnnotation.value,
  selectedSpecificPreset: state.selectedSpecificPreset.value,
  specificRulePresetLoading: Boolean(state.specificRulePresetLoading.value),
  newSpecRule: cloneJsonData(state.newSpecRule),
  newPriorityRule: cloneJsonData(state.newPriorityRule),
  newColorFeat: state.newColorFeat.value,
  newColorVal: state.newColorVal.value,
  newFeatureToAdd: state.newFeatureToAdd.value,
  featureSearch: state.featureSearch.value,
  labelSearch: state.labelSearch.value,
  selectedFeatureIds: Array.from(state.selectedFeatureIds.value || []),
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureSelectionSuppressNextClick: Boolean(
    state.featureSelectionSuppressNextClick.value
  ),
  featureSelectionDrag: cloneJsonData(state.featureSelectionDrag),
  labelReflowLastError: state.labelReflowLastError.value,
  clickedFeature: state.clickedFeature.value,
  clickedPairwiseMatch: state.clickedPairwiseMatch.value,
  clickedLabel: state.clickedLabel.value,
  settingsPopup: state.settingsPopup.value && { ...state.settingsPopup.value },
  featureExtractionPending: Boolean(state.featureExtractionPending.value),
  featureExtractionError: state.featureExtractionError.value,
  featureEditorStatus: cloneJsonData(state.featureEditorStatus),
  matchSequenceSources: cloneJsonData(state.matchSequenceRegistry?.values?.() || []),
  diagramElements: [...state.diagramElements.value],
  diagramElementIds: [...state.diagramElementIds.value],
  diagramElementOriginalTransforms: new Map(
    [...state.diagramElementOriginalTransforms.value].map(([element, transform]) => [
      element,
      cloneJsonData(transform)
    ])
  ),
  legendDragging: Boolean(state.legendDragging.value),
  legendDragStart: cloneJsonData(state.legendDragStart),
  legendOriginalTransform: cloneJsonData(state.legendOriginalTransform.value),
  legendInitialTransform: cloneJsonData(state.legendInitialTransform.value),
  diagramDragging: Boolean(state.diagramDragging.value),
  diagramDragStart: cloneJsonData(state.diagramDragStart),
  lengthBarElement: state.lengthBarElement.value,
  lengthBarOriginalTransform: cloneJsonData(state.lengthBarOriginalTransform.value),
  plotTitleElement: state.plotTitleElement.value,
  plotTitleDragging: Boolean(state.plotTitleDragging.value),
  plotTitleDragStart: cloneJsonData(state.plotTitleDragStart),
  plotTitleAutoTransform: cloneJsonData(state.plotTitleAutoTransform.value)
});

/** @param {ReturnType<typeof captureSessionImportTransientState>} snapshot */
const restoreSessionImportTransientState = (snapshot) => {
  state.suppressCircularMultiRecordDefaults.value =
    snapshot.suppressCircularMultiRecordDefaults;
  state.linearReorderNotice.value = snapshot.linearReorderNotice;
  restoreFoldedLinearSources(state.linearFoldedSources, snapshot.linearFoldedSources);
  restoreRightDrawerState(state, snapshot.rightDrawer);
  state.showCanvasControls.value = snapshot.showCanvasControls;
  state.isPanning.value = snapshot.isPanning;
  Object.assign(state.panStart, cloneJsonData(snapshot.panStart));
  state.selectedAnnotation.value = snapshot.selectedAnnotation;
  state.selectedSpecificPreset.value = snapshot.selectedSpecificPreset;
  state.specificRulePresetLoading.value = snapshot.specificRulePresetLoading;
  replacePlainObject(state.newSpecRule, cloneJsonData(snapshot.newSpecRule));
  replacePlainObject(state.newPriorityRule, cloneJsonData(snapshot.newPriorityRule));
  state.newColorFeat.value = snapshot.newColorFeat;
  state.newColorVal.value = snapshot.newColorVal;
  state.newFeatureToAdd.value = snapshot.newFeatureToAdd;
  state.featureSearch.value = snapshot.featureSearch;
  state.labelSearch.value = snapshot.labelSearch;
  state.selectedFeatureIds.value = new Set(snapshot.selectedFeatureIds);
  state.selectedFeatureAnchorId.value = snapshot.selectedFeatureAnchorId;
  state.featureSelectionStatus.value = snapshot.featureSelectionStatus;
  state.featureSelectionSuppressNextClick.value =
    snapshot.featureSelectionSuppressNextClick;
  Object.assign(
    state.featureSelectionDrag,
    cloneJsonData(snapshot.featureSelectionDrag)
  );
  state.labelReflowLastError.value = snapshot.labelReflowLastError;
  state.clickedFeature.value = snapshot.clickedFeature;
  state.clickedPairwiseMatch.value = snapshot.clickedPairwiseMatch;
  state.clickedLabel.value = snapshot.clickedLabel;
  // The popup of a settings section is open again where it was, as the popups above are.
  if (snapshot.settingsPopup) floatSettingsSection(state, snapshot.settingsPopup.section, snapshot.settingsPopup);
  else returnSettingsSection(state);
  state.featureExtractionPending.value = snapshot.featureExtractionPending;
  state.featureExtractionError.value = snapshot.featureExtractionError;
  Object.assign(state.featureEditorStatus, cloneJsonData(snapshot.featureEditorStatus));
  state.matchSequenceRegistry?.reset?.(cloneJsonData(snapshot.matchSequenceSources));
  state.diagramElements.value = [...snapshot.diagramElements];
  state.diagramElementIds.value = [...snapshot.diagramElementIds];
  state.diagramElementOriginalTransforms.value = new Map(
    snapshot.diagramElementOriginalTransforms
  );
  state.legendDragging.value = snapshot.legendDragging;
  Object.assign(state.legendDragStart, cloneJsonData(snapshot.legendDragStart));
  state.legendOriginalTransform.value = cloneJsonData(snapshot.legendOriginalTransform);
  state.legendInitialTransform.value = cloneJsonData(snapshot.legendInitialTransform);
  state.diagramDragging.value = snapshot.diagramDragging;
  Object.assign(state.diagramDragStart, cloneJsonData(snapshot.diagramDragStart));
  state.lengthBarElement.value = snapshot.lengthBarElement;
  state.lengthBarOriginalTransform.value = cloneJsonData(
    snapshot.lengthBarOriginalTransform
  );
  state.plotTitleElement.value = snapshot.plotTitleElement;
  state.plotTitleDragging.value = snapshot.plotTitleDragging;
  Object.assign(state.plotTitleDragStart, cloneJsonData(snapshot.plotTitleDragStart));
  state.plotTitleAutoTransform.value = cloneJsonData(snapshot.plotTitleAutoTransform);
  state.semanticFileWatchersSuppressed.value =
    snapshot.semanticFileWatchersSuppressed;
  state.skipCaptureBaseConfig.value = snapshot.skipCaptureBaseConfig;
};

// Session Load's rollback holds both drawings; the shown one is restored
// last, so the shown mode's artifact values are its own.
/** @param {DrawingState} drawing */
const captureDrawingImportSnapshot = (drawing) => ({
  config: cloneJsonData(buildConfigData(drawing)),
  ui: cloneJsonData(buildUiStateData(drawing)),
  features: buildFeatureStateData(drawing),
  editorState: buildEditorStateData(drawing),
  orthogroupState: buildOrthogroupStateData(drawing),
  importedComparisonIntent: cloneJsonData(drawing.importedComparisonIntent),
  fileLegendCaptions: Array.from(drawing.fileLegendCaptions.value || [])
});

const captureSessionImportSnapshot = () => ({
  drawings: {
    circular: captureDrawingImportSnapshot(state.drawings.circular),
    linear: captureDrawingImportSnapshot(state.drawings.linear)
  },
  // The cast: outside a Load the shown mode is a drawing's ('circular' or 'linear'), and a Load replaces
  // the saved mode it writes first before it returns (no `await` between).
  ui: cloneJsonData(buildUiStateData(state.drawings[/** @type {DiagramMode} */ (state.mode.value)])),
  richFeaturePopup: Boolean(state.richFeaturePopup.value),
  files: cloneLiveFileState(state.drawings.linear),
  results: state.results.value,
  collinearGroups: state.collinearGroups.value,
  runState: buildRunStateData(),
  losatCache: state.losatCache.value,
  losatDerivedCache: state.losatDerivedCache.value,
  proteinIdentityManifest: toRaw(state.proteinIdentityManifest.value)
    === adoptedProteinIdentityManifest
    ? toRaw(state.proteinIdentityManifest.value)
    : cloneJsonData(state.proteinIdentityManifest.value),
  adoptedProteinIdentityManifest,
  legacyProteinRawCandidates: cloneJsonData(state.legacyProteinRawCandidates.value),
  legacyProteinDerivedEvidence: cloneJsonData(state.legacyProteinDerivedEvidence.value),
  losatCacheInfo: cloneJsonData(state.losatCacheInfo.value),
  committedCanonicalSession: isAdoptedCanonicalSession(committedCanonicalSession)
    ? committedCanonicalSession
    : cloneCanonicalSession(committedCanonicalSession),
  activeSessionResourceTable,
  loadedTrackSettings,
  errorLog: state.errorLog.value,
  resultPanelTab: state.resultPanelTab.value,
  transients: captureSessionImportTransientState()
});

/** @param {ReturnType<typeof captureSessionImportSnapshot>} snapshot */
const restoreSessionImportSnapshot = async (snapshot) => {
  const shownMode = snapshot.ui.mode === 'linear' ? 'linear' : 'circular';
  /** @type {Array<'circular' | 'linear'>} */
  const restoreOrder = shownMode === 'linear' ? ['circular', 'linear'] : ['linear', 'circular'];
  state.sessionImportRollbackInProgress.value = true;
  try {
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    state.mode.value = shownMode;
    restoreOrder.forEach((drawingMode) => {
      const drawing = state.drawings[drawingMode];
      const saved = snapshot.drawings[drawingMode];
      applyCapturedConfigData(drawing, saved.config, { resolveTrackPlacements: false });
      applyUiStateData(drawing, saved.ui);
    });
    applyUiStateData(state.drawings[shownMode], snapshot.ui);
    state.richFeaturePopup.value = snapshot.richFeaturePopup;
    restoreLiveFileState(state.drawings.linear, snapshot.files);
    state.losatCache.value = snapshot.losatCache;
    state.losatDerivedCache.value = snapshot.losatDerivedCache;
    adoptedProteinIdentityManifest = snapshot.adoptedProteinIdentityManifest;
    state.proteinIdentityManifest.value = snapshot.proteinIdentityManifest === adoptedProteinIdentityManifest
      ? snapshot.proteinIdentityManifest
      : cloneJsonData(snapshot.proteinIdentityManifest);
    state.legacyProteinRawCandidates.value = cloneJsonData(snapshot.legacyProteinRawCandidates);
    state.legacyProteinDerivedEvidence.value = cloneJsonData(snapshot.legacyProteinDerivedEvidence);
    state.losatCacheInfo.value = cloneJsonData(snapshot.losatCacheInfo);
    loadedTrackSettings = snapshot.loadedTrackSettings;
    setCommittedCanonicalSession(isAdoptedCanonicalSession(snapshot.committedCanonicalSession)
      ? snapshot.committedCanonicalSession
      : cloneCanonicalSession(snapshot.committedCanonicalSession));
    activeSessionResourceTable = snapshot.activeSessionResourceTable;
    state.skipCaptureBaseConfig.value = true;
    applyResultsData(snapshot.results, snapshot.ui);
    restoreOrder.forEach((drawingMode) => {
      const drawing = state.drawings[drawingMode];
      const saved = snapshot.drawings[drawingMode];
      Object.assign(
        drawing.importedComparisonIntent,
        createImportedComparisonIntentState(),
        cloneJsonData(saved.importedComparisonIntent)
      );
      applyFeatureStateData(drawing, saved.features);
      applyOrthogroupStateData(drawing, saved.orthogroupState);
      applyEditorStateData(drawing, saved.editorState, { normalized: true });
      drawing.fileLegendCaptions.value = new Set(saved.fileLegendCaptions);
    });
    state.collinearGroups.value = snapshot.collinearGroups;
    applyRunStateData(snapshot.runState);
    state.errorLog.value = snapshot.errorLog;
    state.resultPanelTab.value = snapshot.resultPanelTab;
    await nextTick();
    recordSessionLifecycleEvent('session-rollback-source-restored');
    restoreSessionImportTransientState(snapshot.transients);
    recordSessionLifecycleEvent('session-rollback-transients-reconciled');
    await nextTick();
  } finally {
    state.sessionImportRollbackInProgress.value = false;
  }
};

/** @param {Record<string, unknown>} target */
const clearObject = (target) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
};

// Both drawings return to their defaults before a Session installs its own.
const resetSessionBaseline = () => {
  activePreviewRuntime?.clearActiveRuntime?.();
  preservedCliOptions = null;
  setCommittedCanonicalSession(null);
  loadedTrackSettings = null;
  activeSessionResourceTable = null;
  adoptedProteinIdentityManifest = null;
  state.sessionResourceDiscoveryDeferred.value = false;
  resetSettingsState(state);
  resetLayoutState(state);
  resetRightDrawerState(state);
  // Reset Settings keeps which records are ON (PD-OI-070); a new Session draws
  // every record its slices do not turn OFF (PD-OI-091).
  for (const drawing of Object.values(state.drawings)) drawing.recordsOff.splice(0);
  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.sessionTitle.value = '';
  state.errorLog.value = null;
  state.results.value = [];
  state.failedGeneratePreservedResult.value = false;
  state.selectedResultIndex.value = 0;
  state.resultPanelTab.value = 'preview';
  state.lastRunInfo.value = null;
  state.trackSlotResolvedGeometry.value = null;
  state.annotationWarnings.value = [];
  state.featureIdentityNotices.value = [];
  state.featureEditRemovalCount.value = 0;
  state.comparisonWarnings.value = [];
  state.labelFitReports.value = [];
  applyFiles(null, { targetDrawing: state.drawings.linear });
  state.losatCache.value = new Map();
  state.losatDerivedCache.value = new Map();
  state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
  state.legacyProteinRawCandidates.value = { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [] };
  state.legacyProteinDerivedEvidence.value = { schema: LEGACY_LOSAT_DERIVED_CACHE_SCHEMA, entries: [] };
  state.losatCacheInfo.value = [];
  state.orthogroups.value = [];
  state.collinearGroups.value = [];
  state.featureOrthogroupIndex.value = new Map();
  state.selectedOrthogroupId.value = '';
  state.selectedOrthogroupAlignmentFeature.value = '';
  state.extractedFeatures.value = [];
  if (state.biologicalFeatures) state.biologicalFeatures.value = [];
  state.featureRecordIds.value = [];
  state.selectedFeatureRecordIdx.value = 0;
  Object.values(state.drawings).forEach((/** @type {DrawingState} */ drawing) => {
    Object.assign(drawing.importedComparisonIntent, createImportedComparisonIntentState());
    clearObject(drawing.orthogroupNameOverrides);
    clearObject(drawing.orthogroupDescriptionOverrides);
    clearObject(drawing.orthogroupDormantOverrides);
    clearObject(drawing.featureColorOverrides);
    drawing.featureVisibilityManualRules.splice(0);
    clearObject(drawing.featureOverrides);
    clearObject(drawing.featureStrokeOverrides);
    drawing.canonicalLabelOverrideRows.value = [];
    clearObject(drawing.labelTextBulkOverrides);
    drawing.legendEntries.value = [];
    drawing.dormantLegendEntries.value = [];
  });
  state.generatedMode.value = 'circular';
  state.generatedLegendPosition.value = 'left';
  state.generatedMultiRecordCanvas.value = false;
  state.generatedCircularPlotTitlePosition.value = 'none';
};

/**
 * The preview zoom and pan as a Session's `ui` holds them.
 * @typedef {{ zoom: number, canvasPan: { x: number, y: number } }} PersistedPreviewViewport
 */
// The view to save is not always the shown one: Review alignment shows the figure fitted and
// the pan/zoom owner holds the view from before it. The composition root registers that owner's
// reader of the held view; the Session writer and the UI capture with navigation save the held
// view while there is one, and the shown view otherwise.
/** @type {(() => PersistedPreviewViewport | null) | null} */
let heldPreviewViewportReader = null;
/** @param {(() => PersistedPreviewViewport | null) | null} reader */
export const setHeldPreviewViewportReader = (reader) => { heldPreviewViewportReader = reader; };
/** @returns {PersistedPreviewViewport} */
const previewViewportToPersist = () => {
  const { zoom, canvasPan } = heldPreviewViewportReader?.() ?? { zoom: state.zoom.value, canvasPan: state.canvasPan };
  return { zoom, canvasPan: { x: canvasPan.x, y: canvasPan.y } };
};

/** @param {DrawingState} drawing */
export const buildUiStateData = (drawing, { includePreviewNavigation = true } = {}) => {
  const ui = {
    title: String(state.sessionTitle.value || ''),
    mode: state.mode.value,
    canvasPadding: { ...drawing.canvasPadding },
    selectedResultIndex: state.selectedResultIndex.value,
    generatedLegendPosition: state.generatedLegendPosition.value,
    generatedMode: state.generatedMode.value,
    generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
    generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
      state.generatedCircularPlotTitlePosition.value
    ),
    layoutPreferences: cloneJsonData(drawing.layoutPreferences),
    featurePanelTab: state.featurePanelTab.value,
    cInputType: state.cInputType.value,
    lInputType: state.lInputType.value,
    losatProgram: drawing.losatProgram.value,
    downloadDpi: state.downloadDpi.value,
    autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
    linearTypographyLinked: Boolean(drawing.linearTypographyLinked.value),
    paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
    // App-level settings (History and the Load rollback restore them with the UI).
    // The Rich Feature Popup is a UI preference outside History (D-03): the
    // Session and the Load rollback save it on their own.
    losatExecution: buildLosatExecutionData(),
    appliedPaletteName: state.appliedPaletteName.value,
    appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
    pendingPaletteName: drawing.pendingPaletteName.value,
    pendingPaletteColors: cloneColors(drawing.pendingPaletteColors.value),
    legendCurrentOffset: { ...state.legendCurrentOffset },
    diagramOffset: { ...state.diagramOffset },
    lengthBarUserOffset: { ...state.lengthBarUserOffset },
    plotTitleUserOffset: { ...state.plotTitleUserOffset }
  };

  if (includePreviewNavigation) Object.assign(ui, previewViewportToPersist());

  return ui;
};

// `mode` is not restored here: History restores run the mode transition
// first, and Session rollback sets it with the restored document.
/**
 * @param {DrawingState} drawing
 * @param {SavedUiFields} [ui] What `buildUiStateData` wrote (History, the Load rollback).
 * @param {{ restorePreviewNavigation?: boolean }} [options]
 */
export const applyUiStateData = (drawing, ui = {}, { restorePreviewNavigation = true } = {}) => {
  if (typeof ui.title === 'string') state.sessionTitle.value = ui.title;
  if (ui.cInputType) state.cInputType.value = ui.cInputType;
  if (ui.lInputType) state.lInputType.value = ui.lInputType;
  if (ui.losatProgram) {
    const program = String(ui.losatProgram);
    drawing.losatProgram.value = ['blastn', 'tblastx', 'blastp'].includes(program) ? program : 'blastn';
  }
  if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
  state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
  reconcileImportedLinearTypographyLink({
    adv: drawing.adv,
    linked: drawing.linearTypographyLinked,
    ui
  });
  state.paletteInstantPreviewEnabled.value = Boolean(ui.paletteInstantPreviewEnabled);
  // The cast: a plain object, and every field of the type is an optional unknown.
  if (isPlainObject(ui.losatExecution)) applyLosatExecutionData(/** @type {SessionLosatExecution} */ (ui.losatExecution));
  if (ui.featurePanelTab === 'labels' || ui.featurePanelTab === 'colors') {
    state.featurePanelTab.value = ui.featurePanelTab;
  }

  if (ui.generatedMode) state.generatedMode.value = ui.generatedMode === 'linear' ? 'linear' : 'circular';
  if (ui.generatedLegendPosition) {
    state.generatedLegendPosition.value = legendPositionOrDefault(
      ui.generatedLegendPosition, state.generatedMode.value
    );
  }
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    state.generatedMultiRecordCanvas.value = Boolean(ui.generatedMultiRecordCanvas);
  }
  if (ui.generatedCircularPlotTitlePosition || ui.circularPlotTitlePosition) {
    state.generatedCircularPlotTitlePosition.value = hasStoredLayoutValue(ui.generatedCircularPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition)
      : normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition);
  }

  restorePaletteStateFromSession(drawing, ui);
  restoreLayoutPreferences(drawing, ui);

  if (ui.legendCurrentOffset) {
    state.legendCurrentOffset.x = Number(ui.legendCurrentOffset.x) || 0;
    state.legendCurrentOffset.y = Number(ui.legendCurrentOffset.y) || 0;
  }
  if (ui.diagramOffset) {
    state.diagramOffset.x = Number(ui.diagramOffset.x) || 0;
    state.diagramOffset.y = Number(ui.diagramOffset.y) || 0;
  }
  if (ui.lengthBarUserOffset) {
    state.lengthBarUserOffset.x = Number(ui.lengthBarUserOffset.x) || 0;
    state.lengthBarUserOffset.y = Number(ui.lengthBarUserOffset.y) || 0;
  }
  if (ui.plotTitleUserOffset) {
    state.plotTitleUserOffset.x = Number(ui.plotTitleUserOffset.x) || 0;
    state.plotTitleUserOffset.y = Number(ui.plotTitleUserOffset.y) || 0;
  }

  if (ui.canvasPadding) Object.assign(drawing.canvasPadding, readCanvasPadding(ui.canvasPadding).value);
  if (restorePreviewNavigation) {
    if (ui.canvasPan) {
      state.canvasPan.x = Number(ui.canvasPan.x) || 0;
      state.canvasPan.y = Number(ui.canvasPan.y) || 0;
    }
    if (typeof ui.zoom === 'number') state.zoom.value = ui.zoom;
  }

};

/**
 * @param {unknown} [resultsData] Results: committed ones, or imported name/content pairs.
 * @param {{ selectedResultIndex?: unknown }} [ui]
 */
export const applyResultsData = (resultsData = [], ui = {}) => {
  state.failedGeneratePreservedResult.value = false;
  if (Array.isArray(resultsData)) {
    const logicalResults = resultsData.every(isCommittedSvgResult)
      ? resultsData
      : normalizeLogicalResults(resultsData.map((res, idx) => (
          isCommittedSvgResult(res)
            ? res
            : { name: res?.name || `Result ${idx + 1}`, content: res?.content || '' }
        )));
    const committedCount = logicalResults.filter(isCommittedSvgResult).length;
    if (committedCount !== 0 && committedCount !== logicalResults.length) {
      throw new Error('Committed and imported SVG Results cannot be mixed.');
    }
    // Results restored without their committed metadata (a History
    // checkpoint) carry the feature types of the committed request (R13).
    state.results.value = committedCount === logicalResults.length
      ? logicalResults
      : admitLegacyImportedResults(createLegacyImportResultSource(logicalResults), {
          selectedFeatureTypes: committedCanonicalSession?.renderRequest?.diagramOptions?.selectedFeaturesSet
        });
  } else {
    state.results.value = [];
  }

  const resultCount = state.results.value.length;
  if (resultCount > 0) {
    const stored = ui.selectedResultIndex;
    // `typeof` adds nothing to `Number.isInteger`; it narrows `stored`.
    const desiredIndex = typeof stored === 'number' && Number.isInteger(stored) && stored >= 0 ? stored : 0;
    state.selectedResultIndex.value = Math.min(desiredIndex, resultCount - 1);
  } else {
    state.selectedResultIndex.value = 0;
  }
};

// The catalog features and groups are read from their raw arrays: the copy
// (Session, History checkpoint, import rollback) holds the same values, without
// a reactive read of each of their 10k-40k nested properties.
/** @param {DrawingState} drawing */
export const buildFeatureStateData = (drawing) => ({
  extractedFeatures: sanitizeExtractedFeaturesForSession(toRaw(state.extractedFeatures.value)),
  biologicalFeatures: sanitizeExtractedFeaturesForSession(toRaw(state.biologicalFeatures?.value)),
  featureRecordIds: cloneJsonData(state.featureRecordIds.value),
  selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
  featureColorOverrides: cloneJsonData(drawing.featureColorOverrides),
  featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(drawing.featureVisibilityManualRules),
  featureOverrides: featureOverridesForState(drawing.featureOverrides),
  labelOverrideRows: cloneJsonData(drawing.canonicalLabelOverrideRows.value),
  labelTextBulkOverrides: cloneJsonData(drawing.labelTextBulkOverrides)
});

// The catalog features of the shown Result (artifact values).
/** @param {SavedFeatureFields} [features] */
const applyFeatureArtifactData = (features = {}) => {
  // The casts: `Array.isArray` proves an array, not its elements.
  state.extractedFeatures.value = Array.isArray(features.extractedFeatures)
    ? /** @type {unknown[]} */ (features.extractedFeatures)
    : [];
  if (state.biologicalFeatures) {
    state.biologicalFeatures.value = Array.isArray(features.biologicalFeatures)
      ? /** @type {unknown[]} */ (features.biologicalFeatures)
      : [];
  }
  state.featureRecordIds.value = Array.isArray(features.featureRecordIds)
    ? /** @type {unknown[]} */ (features.featureRecordIds)
    : [];
};

/**
 * @param {DrawingState} drawing
 * @param {SavedFeatureFields} [features]
 */
export const applyFeatureStateData = (drawing, features = {}) => {
  applyFeatureArtifactData(features);
  state.selectedFeatureRecordIdx.value = Number.isInteger(features.selectedFeatureRecordIdx)
    ? features.selectedFeatureRecordIdx
    : 0;
  applyDrawingFeatureData(drawing, features);
};

// A drawing's per-feature edits as saved (a Session 46 slice's `features`, a
// flat draft's, or a History checkpoint's).
/**
 * @param {DrawingState} drawing
 * @param {SavedFeatureFields} [features]
 */
const applyDrawingFeatureData = (drawing, features = {}) => {
  // The readers Load lists what they cannot restore with (`sliceFieldRepairs`).
  replacePlainObject(drawing.featureColorOverrides, cloneJsonObject(readFeatureColorEdits(features.featureColorOverrides).value));
  drawing.featureVisibilityManualRules.splice(
    0,
    drawing.featureVisibilityManualRules.length,
    ...normalizeFeatureVisibilityRulesForSession(readList(features.featureVisibilityManualRules).value)
  );
  replacePlainObject(drawing.featureOverrides, featureOverridesForState(features.featureOverrides));
  drawing.canonicalLabelOverrideRows.value = cloneJsonData(readLabelOverrideRows(features.labelOverrideRows).value);
  replacePlainObject(drawing.labelTextBulkOverrides, readStringMap(features.labelTextBulkOverrides).value);
};

/** @param {DrawingState} drawing */
export const buildOrthogroupStateData = (drawing) => ({
  groups: Array.isArray(state.orthogroups.value) ? cloneJsonData(toRaw(state.orthogroups.value)) : [],
  selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
  orthogroupNameOverrides: cloneStringMap(drawing.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: cloneStringMap(drawing.orthogroupDescriptionOverrides),
  orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(drawing.orthogroupDormantOverrides)
});

export const buildRunStateData = () => ({
  lastRunInfo: cloneJsonData(state.lastRunInfo.value),
  pairwiseMatchFactors: cloneJsonObject(state.pairwiseMatchFactors.value)
});

/** @param {NonNullable<UnvalidatedFields<'lastRunInfo' | 'pairwiseMatchFactors'>>} [runState] */
export const applyRunStateData = (runState = {}) => {
  state.lastRunInfo.value = runState.lastRunInfo ? cloneJsonData(runState.lastRunInfo) : null;
  state.pairwiseMatchFactors.value = cloneJsonObject(runState.pairwiseMatchFactors);
};

/**
 * @param {{ status?: string, generationId?: string, error?: string | null, summaryCount?: number,
 *   detailsCacheSize?: number }} [updates]
 */
const setFeatureEditorStatusData = (updates = {}) => {
  if (!state.featureEditorStatus || typeof state.featureEditorStatus !== 'object') return;
  Object.assign(state.featureEditorStatus, {
    status: updates.status ?? state.featureEditorStatus.status,
    generationId: updates.generationId ?? state.featureEditorStatus.generationId,
    error: updates.error === undefined ? state.featureEditorStatus.error : updates.error,
    summaryCount: updates.summaryCount ?? state.featureEditorStatus.summaryCount,
    detailsCacheSize: updates.detailsCacheSize ?? state.featureEditorStatus.detailsCacheSize
  });
};

const synchronizeRestoredFeatureSummaryStatus = ({ generationId = 'session-load' } = {}) => {
  const summaryCount = Array.isArray(state.extractedFeatures.value)
    ? state.extractedFeatures.value.length
    : 0;
  if (summaryCount === 0) return false;
  state.featureExtractionPending.value = false;
  state.featureExtractionError.value = null;
  setFeatureEditorStatusData({
    status: 'summary-ready',
    generationId,
    error: null,
    summaryCount
  });
  return true;
};

/**
 * @param {DrawingState} drawing
 * @param {{ status?: string, recoveredFeatureState?: SavedFeatureFields | null,
 *   migratedEditorState?: Record<string, unknown> | null, warning?: string | null } | null | undefined} plan
 *   What this reads of a feature recovery plan (`buildSessionFeatureRecoveryPlan`, or the failed plan Load builds).
 * @param {{ generationId?: string }} [options]
 */
const applySessionFeatureRecoveryPlan = (drawing, plan, { generationId = 'session-feature-recovery' } = {}) => {
  state.featureExtractionPending.value = false;

  if (plan?.status === 'recovered' || plan?.status === 'aligned') {
    if (plan.recoveredFeatureState) applyFeatureStateData(drawing, plan.recoveredFeatureState);
    if (plan.migratedEditorState) applyEditorStateData(drawing, plan.migratedEditorState);
    state.featureExtractionError.value = null;
    setFeatureEditorStatusData({
      status: 'summary-ready',
      generationId,
      error: plan.warning ? tm(plan.warning) : null,
      // The cast: `Array.isArray` on the line above; the elements stay unvalidated.
      summaryCount: Array.isArray(plan.recoveredFeatureState?.extractedFeatures)
        ? /** @type {unknown[]} */ (plan.recoveredFeatureState.extractedFeatures).length
        : state.extractedFeatures.value.length
    });
    return;
  }

  if (plan?.status === 'unrecoverable' || plan?.status === 'failed') {
    const warning = plan.warning ? tm(plan.warning) : t('Feature metadata recovery failed. The SVG preview remains available.');
    state.featureExtractionError.value = { summary: warning, details: [] };
    setFeatureEditorStatusData({
      status: 'failed',
      generationId,
      error: warning,
      summaryCount: 0
    });
    return;
  }

  synchronizeRestoredFeatureSummaryStatus({ generationId });
};

// E1: one mode's generated artifact as Save reads it: the installed artifact
// of the displayed mode, or the other mode's slot, which the composition root
// keeps while that mode is not shown.
/**
 * @typedef {Readonly<Record<string, any>> & {
 *   mode: 'circular' | 'linear',
 *   committedCanonicalSession: any,
 *   activeSessionResourceTable: any
 * }} SessionArtifactView
 */
/** @returns {SessionArtifactView} */
const displayedArtifactView = () => Object.freeze({
  ...Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [key, state[key].value])),
  originalLegendOrder: state.originalLegendOrder.value,
  mode: state.generatedMode.value === 'linear' ? 'linear' : 'circular',
  committedCanonicalSession,
  activeSessionResourceTable
});
/**
 * @param {Readonly<ArtifactSlot>} slot
 * @returns {SessionArtifactView}
 */
const stashedArtifactView = (slot) => Object.freeze({
  ...slot.values,
  originalLegendOrder: slot.legendInventory,
  mode: slot.mode,
  committedCanonicalSession: slot.runtimeState?.canonical?.committedCanonicalSession ?? null,
  activeSessionResourceTable: slot.runtimeState?.canonical?.activeSessionResourceTable ?? null
});
/** @param {SessionArtifactView | null} view */
const artifactHasResult = (view) => Array.isArray(view?.results) && view.results.length > 0;
// Save writes every Result (Owner, 2026-10-07). The top level holds the shown
// mode's set when it has a Result, otherwise the other mode's set;
// `otherModeResult` holds the remaining set only when it has a Result.
/** @param {Readonly<ArtifactSlot> | null} otherSlot */
const chooseSessionArtifacts = (otherSlot) => {
  const displayed = displayedArtifactView();
  const other = otherSlot && otherSlot.mode !== displayed.mode ? stashedArtifactView(otherSlot) : null;
  if (!other || !artifactHasResult(other)) return { artifact: displayed, otherArtifact: null };
  if (!artifactHasResult(displayed)) return { artifact: other, otherArtifact: null };
  return { artifact: displayed, otherArtifact: other };
};
// Live edits are committed into a Result's content at once, so its content is
// what Save writes, also for the displayed Result.
/** @param {SessionArtifactView} view */
const serializeArtifactResults = (view) => normalizeLogicalResults(view.results.map(
  (/** @type {Record<string, any>} */ res, /** @type {number} */ idx) => ({
    name: res.name || `Result ${idx + 1}`,
    content: res.content
  })
));
// The artifact's own catalog, generated Legend inventory and colors, stroke
// defaults and alignment Reset receipt; the Legend and stroke edits are its
// mode's drawing's (Session 46 `modes.<m>.editorState`).
/** @param {SessionArtifactView} view */
const buildArtifactEditorState = (view) => {
  return {
    legend: {
      originalOrder: cloneJsonArray(view.originalLegendOrder),
      originalColors: cloneStringMap(view.originalLegendColors)
    },
    originalSvgStroke: {
      color: view.originalSvgStroke?.color ?? null,
      width: view.originalSvgStroke?.width ?? null
    },
    alignmentResetReceipt: cloneJsonValue(view.similarityAlignmentResetReceipt, null),
    featureCatalog: admittedFeatureCatalog(view.featureCatalog)
  };
};
// D-25 (PD-OI-079): Results need their admitted catalog; a Result without one
// (a legacy Session) is saved after a Generate in its mode. The error names that
// mode when it is not the one shown, and its Generate action runs there.
/**
 * @param {FeatureCatalog | null} featureCatalog The artifact's admitted catalog.
 * @param {FeatureCatalogResult[]} logicalResults
 * @param {'circular' | 'linear'} mode
 * @returns {FeatureCatalog | null} The catalog checked against the Results; null without Results.
 */
const admitSavedFeatureCatalog = (featureCatalog, logicalResults, mode) => {
  if (logicalResults.length === 0) return null;
  const requiresGenerate = () => sessionSaveRequiresGenerate(
    mode === state.mode.value ? {} : { diagramMode: mode }
  );
  if (!featureCatalog) throw requiresGenerate();
  try {
    return validateFeatureCatalog(
      featureCatalog,
      logicalResults,
      { adopt: true, mode }
    );
  } catch (error) {
    console.warn('Session feature catalog validation failed.', normalizeUserFacingError(error));
    throw requiresGenerate();
  }
};
/** @param {SessionArtifactView} view */
const artifactRunMetadata = (view) => ({
  ...(view.trackSlotResolvedGeometry
    ? { trackSlotGeometry: cloneJsonData(view.trackSlotResolvedGeometry) } : {}),
  annotationWarnings: cloneJsonData(view.annotationWarnings),
  ...(view.featureIdentityNotices?.length
    ? { featureIdentityNotices: cloneJsonData(view.featureIdentityNotices) } : {}),
  ...(view.comparisonWarnings?.length
    ? { comparisonWarnings: cloneJsonData(view.comparisonWarnings) } : {})
});
/** @param {SessionArtifactView} view */
const artifactCliInvocation = (view) => {
  const invocation = view.lastRunInfo?.invocation;
  return isCliInvocationSessionExportable(invocation) ? cloneJsonData(invocation) : undefined;
};
// The generated layout and palette of an artifact's Results (`ui` fields).
/** @param {SessionArtifactView} view */
const artifactUiState = (view) => ({
  selectedResultIndex: view.selectedResultIndex,
  generatedLegendPosition: view.generatedLegendPosition,
  generatedMultiRecordCanvas: Boolean(view.generatedMultiRecordCanvas),
  generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
    view.generatedCircularPlotTitlePosition
  ),
  appliedPaletteName: view.appliedPaletteName,
  appliedPaletteColors: cloneColors(view.appliedPaletteColors)
});
// The other mode's set: its Results, admitted catalog and committed request
// are written as one unit beside the top-level set (PD-OI-045).
/** @param {SessionArtifactView} view */
const prepareOtherModeArtifact = (view) => {
  const results = serializeArtifactResults(view);
  const editorState = buildArtifactEditorState(view);
  editorState.featureCatalog = admitSavedFeatureCatalog(editorState.featureCatalog, results, view.mode);
  if (!view.committedCanonicalSession) throw sessionSaveRequiresGenerate({ diagramMode: view.mode });
  const committed = isAdoptedCanonicalSession(view.committedCanonicalSession)
    ? view.committedCanonicalSession
    : cloneCanonicalSession(view.committedCanonicalSession);
  return { view, results, editorState, committed };
};
/**
 * @param {ReturnType<typeof prepareOtherModeArtifact>} other
 * @param {CanonicalRenderRequest} renderRequest Its resource references name the Session's table.
 * @returns {SessionOtherModeResult}
 */
const buildOtherModeResult = (other, renderRequest) => {
  const cliInvocation = artifactCliInvocation(other.view);
  return {
    renderRequest,
    results: other.results,
    editorState: {
      featureCatalog: other.editorState.featureCatalog,
      alignmentResetReceipt: other.editorState.alignmentResetReceipt ?? null,
      legend: {
        originalOrder: other.editorState.legend.originalOrder,
        originalColors: other.editorState.legend.originalColors
      },
      originalSvgStroke: other.editorState.originalSvgStroke
    },
    ui: artifactUiState(other.view),
    runMetadata: artifactRunMetadata(other.view),
    ...(cliInvocation ? { cliInvocation } : {})
  };
};

// `askChoice`: the root's choice dialog; Save asks with it before a repeated or a large download.
/**
 * @typedef {{
 *   askChoice?: DownloadChoicePort,
 *   drawing: DrawingState,
 *   linearRecordCatalog?: any,
 *   omittedAnnotationRecordKeys?: string[],
 *   recordDisplayRows?: any,
 *   modes: Record<'circular' | 'linear', SessionModeSlice & { config: ActiveWebConfig }>,
 *   savedUi: Record<string, any>,
 *   isCurrent: () => boolean,
 *   artifact: SessionArtifactView,
 *   otherArtifact: SessionArtifactView | null
 * }} ExportSessionDocumentOptions
 */

/**
 * How a Save ended: the download was handed off, or the person or a newer
 * operation canceled it.
 * @typedef {{ status: 'saved', blob: Blob, filename: string }
 *   | { status: 'canceled', compressedSize?: number }} SessionSaveResult
 */

/**
 * @param {string | null | undefined} titleOverride
 * @param {ExportSessionDocumentOptions} options
 * @returns {Promise<SessionSaveResult>}
 */
const exportSessionDocument = async (
  titleOverride = null,
  {
    askChoice = async () => 'cancel',
    drawing, linearRecordCatalog = null, omittedAnnotationRecordKeys = [], recordDisplayRows = null, modes, savedUi,
    isCurrent, artifact, otherArtifact
  }
) => {
  const resolvedTitle =
    typeof titleOverride === 'string'
      ? titleOverride.trim()
      : typeof state.sessionTitle?.value === 'string'
        ? state.sessionTitle.value.trim()
        : '';
  const sessionFilename = buildSessionFilename(resolvedTitle);
  if (lastSessionFilename && lastSessionFilename === sessionFilename) {
    const choice = await askChoice({
      title: t('Download the Session again?'),
      description: t('"{file}" was already downloaded. Your browser may replace that file or save this one under a new name.', { file: sessionFilename }),
      actions: [{ id: 'download', label: t('Download again'), tone: 'primary' }]
    });
    if (choice !== 'download' || !isCurrent()) {
      recordSessionLifecycleEvent('session-save-download-canceled', {
        reason: 'repeat-download'
      });
      return { status: 'canceled' };
    }
  }

  recordSessionLifecycleEvent('session-save-projection-start');
  const logicalResults = serializeArtifactResults(artifact);
  const editorState = buildArtifactEditorState(artifact);
  editorState.featureCatalog = admitSavedFeatureCatalog(editorState.featureCatalog, logicalResults, artifact.mode);
  const other = otherArtifact ? prepareOtherModeArtifact(otherArtifact) : null;

  const {
    entries: losatEntries,
    validatedManifest,
    manifestValidated
  } = serializeLosatCache();
  const exportableCliInvocation = artifactCliInvocation(artifact);
  // Each slice's unmanaged overrides are checked against its own mode (PD-OI-009, OV-106).
  for (const mode of SLICE_MODES) {
    const config = modes[mode].config;
    Object.assign(config.adv, normalizedArrowGeometryState(config.adv));
    config.unmanagedConfigOverrides = await validateUnmanagedConfigOverrides({
      mode,
      configOverrides: config.unmanagedConfigOverrides,
      source: 'draft'
    });
  }
  let committed = isAdoptedCanonicalSession(artifact.committedCanonicalSession)
    ? artifact.committedCanonicalSession
    : cloneCanonicalSession(artifact.committedCanonicalSession);
  const settingsOnly = !artifact.committedCanonicalSession && logicalResults.length === 0
    && !hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs });
  if (committed) {
    try {
      const adoptedCommitted = isAdoptedCanonicalSession(committed);
      const projected = projectCanonicalSessionRequest({
        ...committed,
        sessionResourceTable: adoptedCommitted ? artifact.activeSessionResourceTable : null,
        deferResourceContent: adoptedCommitted,
        adoptCanonicalPayloads: adoptedCommitted
      });
      validateCurrentWriterActiveConfig({
        mode: projected.mode,
        storedConfig: modes[projected.mode === 'linear' ? 'linear' : 'circular'].config
      });
    } catch (error) {
      console.warn('Session active configuration validation failed.', normalizeUserFacingError(error));
      throw recognizedCauseOr(error, diagnosticError('INPUT_INVALID', { field: 'config', reason: 'FIELDS' }));
    }
  }
  if (!committed && !settingsOnly) {
    const comparisonPlanSnapshot = state.mode.value === 'linear'
      ? resolveLinearComparisonPlan({
          plan: drawing.linearComparisonPlan,
          sequences: normalizeLinearSeqList(state.linearSeqs),
          recordsOff: drawing.recordsOff,
          layout: drawing.linearRecordLayoutEnabled?.value
            ? drawing.linearRecordRows
            : [],
          losatProgram: drawing.losatProgram?.value,
          blastpMode: drawing.losat?.blastp?.mode
        })
      : null;
    const activeFiles = await serializeActiveRenderFiles(
      state.mode.value,
      state,
      drawing,
      {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog
      }
    );
    committed = buildCanonicalRenderRequest({
      state,
      drawing,
      filesData: activeFiles,
      recordDisplayRows: recordDisplayRows?.value || [],
      comparisonPlanSnapshot,
      omittedAnnotationRecordKeys,
      omitUnsourcedDepth: true
    });
  }
  const canonical = await assembleSessionResources(state, committed, state.drawings.linear, other?.committed ?? null);
  await validateSimilarityAlignmentResetReceipt(editorState.alignmentResetReceipt, canonical);
  if (other) {
    await validateSimilarityAlignmentResetReceipt(other.editorState.alignmentResetReceipt, {
      renderRequest: canonical.otherRenderRequest, resources: canonical.resources
    });
  }
  const legacyRawCandidates = serializableLegacyProteinCandidateEnvelope(
    state.legacyProteinRawCandidates.value
  );
  const legacyDerivedEvidence = normalizeLegacyDerivedEvidence(
    state.legacyProteinDerivedEvidence.value
  );
  if (
    (!manifestValidated || validatedManifest !== state.proteinIdentityManifest.value)
    && !validateProteinIdentityManifest(state.proteinIdentityManifest.value)
  ) {
    throw new Error('Save Session requires a valid protein identity manifest.');
  }
  /** @type {GbdrawSession} */
  const sessionData = /** @satisfies {GbdrawSession} */ ({
    format: 'gbdraw-session',
    version: SESSION_VERSION,
    createdAt: new Date().toISOString(),
    title: resolvedTitle || undefined,
    ui: savedUi,
    renderRequest: canonical.renderRequest,
    resources: canonical.resources,
    webFiles: canonical.webFiles,
    results: logicalResults,
    ...(!settingsOnly ? { runMetadata: artifactRunMetadata(artifact) } : {}),
    ...(other ? { otherModeResult: buildOtherModeResult(other, canonical.otherRenderRequest) } : {}),
    editorState,
    // The group names are the Linear drawing's (`modes.linear.config.webEdits`).
    orthogroupState: {
      selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
      selectedOrthogroupAlignmentFeature: String(state.selectedOrthogroupAlignmentFeature.value || '')
    },
    modes,
    ...(preservedCliOptions ? { cliOptions: cloneJsonData(preservedCliOptions) } : {}),
    losatCache: {
      entries: losatEntries
    },
    losatDerivedCache: {
      entries: []
    },
    proteinIdentityManifest: toRaw(state.proteinIdentityManifest.value)
      === adoptedProteinIdentityManifest
      ? toRaw(state.proteinIdentityManifest.value)
      : cloneJsonData(state.proteinIdentityManifest.value),
    cliInvocation: exportableCliInvocation
  });

  const legacyArtifacts = buildSessionLegacyArtifacts({
    legacyRawCandidates,
    legacyDerivedEvidence
  });
  if (legacyArtifacts) {
    sessionData.legacyArtifacts = legacyArtifacts;
  }
  try {
    validateSessionAuthorityInventory(sessionData, SESSION_VERSION);
  } catch (error) {
    console.error('Session writer validation failed.', normalizeUserFacingError(error));
    throw recognizedCauseOr(error, diagnosticError('INPUT_INVALID', { field: 'schema' }));
  }

  recordSessionLifecycleEvent('session-save-projection-end');
  recordSessionLifecycleEvent('session-save-compression-start');
  const compressed = await compressSessionData(sessionData);
  recordSessionLifecycleEvent('session-save-compression-end', {
    compressedSize: compressed.size
  });
  if (!isCurrent()) return { status: 'canceled' };
  const largeMegabytes = largeSessionDownloadMegabytes(compressed);
  if (largeMegabytes !== null && !await askBeforeLargeDownload(askChoice, {
    title: t('Download a large Session?'),
    description: t('The Session file is {size} MB.', { size: largeMegabytes.toFixed(1) })
  })) {
    recordSessionLifecycleEvent('session-save-download-canceled', {
      reason: 'large-download',
      compressedSize: compressed.size
    });
    return { status: 'canceled', compressedSize: compressed.size };
  }
  if (!isCurrent()) return { status: 'canceled' };
  downloadBlob(compressed, sessionFilename);
  recordSessionLifecycleEvent('session-save-download-handoff-completed', {
    compressedSize: compressed.size
  });
  lastSessionFilename = sessionFilename;
  return { status: 'saved', blob: compressed, filename: sessionFilename };
};

/**
 * Asks, before a Session replaces the current one, about the fonts its
 * drawings name (each drawing's `adv.font_family`); false loads nothing.
 * @typedef {(fontFamilies: readonly unknown[]) => Promise<boolean>} SessionFontConfirmation
 */

/**
 * A Result set whose Results Load cannot admit or compose (A20, OV-789):
 * Load reads the file again with that set's Results dropped and listed
 * (`readSessionAuthorityParts`, E-OD16, E-OD17). It never leaves
 * `importSessionDocument`.
 */
class UncomposedResultSet extends Error {
  /** @param {ResultSetPath} at @param {unknown} cause */
  constructor(at, cause) {
    super('RESULT_UNCOMPOSED', { cause });
    this.name = 'UncomposedResultSet';
    this.at = at;
  }
}

/**
 * @param {{ target: { files: ArrayLike<File>, value?: string } }} e
 * @param {Record<string, any>} [options] The options of `importSession` with
 * `isCurrent` and `signal`; `transformLegacyResultSvg` is the composition
 * root's `LegacyResultSvgTransform` port, and `confirmSessionFonts` its
 * `SessionFontConfirmation` port; `uncomposedSets` (Load's own) names the
 * sets a first reading found it cannot compose.
 */
const importSessionDocument = async (e, options = {}) => {
  /** @type {readonly ResultSetPath[]} */
  const uncomposedSets = options.uncomposedSets ?? [];
  const file = e.target.files[0];
  if (!file) return { status: 'skipped' };
  recordSessionLifecycleEvent('sessionSelection');

  const semanticFileWatchersSuppressedBeforeImport = Boolean(
    state.semanticFileWatchersSuppressed.value
  );
  const rollbackStateExtension = options?.rollbackState;
  /** @type {ReturnType<typeof captureSessionImportSnapshot> | null} */
  let rollbackSnapshot = null;
  let rollbackExtensionSnapshot;
  let commitStarted = false;
  const previousAlert = state.errorLog.value;

  try {
    const candidate = await importSessionFile(file, { signal: options.signal });
    if (!options.isCurrent()) return { status: 'canceled' };
    recordSessionLifecycleEvent('session-import-codec-completed', candidate.timings);
    const parsed = candidate.data;
    await assertSafeObjectKeysForImport(parsed);
    if (isLegacyConfigPayload(parsed)) {
      // A legacy configuration names no mode: it is the shown mode's settings.
      // It is read, confirmed and then applied as one step that a failure
      // rolls back (OV-710); its History is not restarted.
      const legacyMode = state.mode.value === 'linear' ? 'linear' : 'circular';
      const legacyConfig = readLegacyConfigPayload(parsed, state.drawings[legacyMode], legacyMode);
      const legacyRepairs = [...legacyConfig.repairs];
      const legacyOverrides = readPlainObject(legacyConfig.value.unmanagedConfigOverrides).value;
      if (Object.keys(legacyOverrides).length > 0) {
        legacyConfig.value.unmanagedConfigOverrides = await readStoredUnmanagedConfigOverrides({
          mode: legacyMode, configOverrides: legacyOverrides, source: 'draft'
        }, 'unmanagedConfigOverrides', legacyRepairs);
      }
      const legacyDialog = describeLoadRepairs(legacyRepairs);
      if (legacyDialog) {
        const confirmLoadRepairs = /** @type {ConfirmLoadRepairs} */ (options.confirmLoadRepairs);
        if (!await confirmLoadRepairs(legacyDialog) || !options.isCurrent()) return { status: 'canceled' };
      }
      rollbackSnapshot = captureSessionImportSnapshot();
      if (typeof rollbackStateExtension?.capture === 'function') {
        rollbackExtensionSnapshot = rollbackStateExtension.capture();
      }
      commitStarted = true;
      applyLegacyConfigPayload(state.drawings[legacyMode], legacyConfig.value);
      return { status: 'legacy', repairs: Object.freeze(legacyRepairs) };
    }

    recordSessionLifecycleEvent('current-session-preflight-start');
    // What Load cannot restore as stored, listed in the dialog before the
    // commit (BA-D25); frozen there, so no step after the commit adds one.
    /** @type {LoadRepair[]} */
    const repairs = [];
    // The one read of the document (D1c 2.2): its band's description, repaired.
    // The document as the file holds it (R4 reads it), and as its band reads it.
    const storedDocument = readSessionInput(parsed);
    const document = readSessionDocument(storedDocument);
    repairs.push(...document.repairs);
    const preflight = await preflightSessionImport(storedDocument, document.value, repairs, uncomposedSets);
    // Preflight read the receipts against the requests Load read (RK-2).
    const { data } = preflight;
    recordSessionLifecycleEvent('current-session-preflight-end');
    const {
      sourceSessionVersion,
      canonicalProjection,
      canonicalProjectionInput,
      restoredConfig,
      projectionResult,
      adoptedCanonicalSession,
      currentResourceTable,
      otherModeCatalog,
      comparisonClassification,
      unmanagedConfigValidation,
      settingsOnly
    } = preflight;
    if (restoredConfig && unmanagedConfigValidation) {
      // A stored draft's overrides are read; a request's are checked.
      restoredConfig.unmanagedConfigOverrides = unmanagedConfigValidation.source === 'draft'
        ? await readStoredUnmanagedConfigOverrides(unmanagedConfigValidation, sourceSessionVersion >= MODE_SCOPED_SESSION_VERSION
          ? `modes.${unmanagedConfigValidation.mode}.config.unmanagedConfigOverrides` : 'config.unmanagedConfigOverrides', repairs)
        : await validateUnmanagedConfigOverrides(unmanagedConfigValidation);
    }
    const canonicalSession = projectionResult !== null;
    const currentSchemaSession = sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION;
    const committedMode = projectionResult?.renderState.mode;
    const savedCurrentWriterMode = sourceSessionVersion >= TYPED_DRAFT_SESSION_VERSION
      ? data.ui?.mode
      : null;
    const ui = canonicalSession
      ? {
          ...projectionResult.editorMetadata.ui,
          ...projectionResult.artifactState.ui,
          mode: ['circular', 'linear'].includes(savedCurrentWriterMode)
            ? savedCurrentWriterMode
            : committedMode
        }
      : (data.ui || {});
    /** @type {'circular' | 'linear'} */
    const candidateMode = readSessionUiMode(ui.mode).value === 'linear' ? 'linear' : 'circular';
    // The Result's mode: its request's. A Session 27-30 has no request; its
    // last Generate's invocation names it, else `ui.mode`, as the CLI reads it
    // (`session_io.session_mode`). 0.13.0 kept a Result across a mode switch,
    // and its `ui.mode` is the shown mode (OV-677). Sessions 29 (tag 0.12.0)
    // and 30 (tag 0.13.0) write `cliInvocation` only for a Generate the CLI
    // can replay; without it the Result reads as the shown mode (PD-OI-096). The flat draft's settings
    // stay with the shown mode (the split's profile mode); its Result edits
    // and the Result's set go with the Result's mode.
    const invocationMode = sourceSessionVersion < CANONICAL_SESSION_MIN_VERSION
      ? /** @type {UnvalidatedFields<'mode'>} */ (data.cliInvocation)?.mode : null;
    /** @type {'circular' | 'linear'} */
    const generatedMode = (committedMode || (invocationMode === 'linear' || invocationMode === 'circular' ? invocationMode : null)
      || ui.mode) === 'linear' ? 'linear' : 'circular';
    // A Session 27-30 saved after a mode switch keeps the shown mode's form (OV-677).
    const formIsResultMode = generatedMode === candidateMode || canonicalSession;
    // The legend position the shown Result was drawn with (OV-652, BA-D11);
    // when the form is not the Result's, it is that Result's side (OV-964).
    const generatedLegend = ui.generatedLegendPosition
      ? readLegendPositionOfMode(ui.generatedLegendPosition, generatedMode, 'ui.generatedLegendPosition',
        formIsResultMode ? 'legend-position-defaulted' : 'legend-side-defaulted') : null;
    if (generatedLegend) repairs.push(...generatedLegend.repairs);
    // PD-OI-086: each mode's drawing gets its own draft. A Session 46 keeps it
    // in its slice; the committed mode's slice lies over its request's
    // projection (`restoredConfig`, plan 4.1). A Session 27-44 has one flat
    // draft, split below once its edits are migrated.
    const modeScopedSession = sourceSessionVersion >= MODE_SCOPED_SESSION_VERSION;
    // A Session 46 committed mode with a Result: its slice over its request's
    // projection (`restoredConfig`); any other slice over its mode's defaults.
    const committedSliceProjected = modeScopedSession && canonicalSession && !settingsOnly && Boolean(restoredConfig);
    // Any other Session 46 slice is stored settings: its track stacks install as saved (R11).
    /** @param {'circular' | 'linear'} mode */
    const isStoredSlice = (mode) => modeScopedSession && !(mode === generatedMode && committedSliceProjected);
    /** @param {'circular' | 'linear'} mode @returns {Record<string, any> | null} */
    const draftConfigOf = (mode) => (modeScopedSession
      ? (mode === generatedMode && committedSliceProjected ? restoredConfig : overlayModeSliceConfig(null, data.modes?.[mode]?.config))
      : restoredConfig);
    const linearDraftConfig = draftConfigOf('linear');
    const circularDraftConfig = draftConfigOf('circular');
    // Each mode's input type (OV-594): the stored one, else its request's
    // (which wins without a row over a stored value of no type), else GenBank.
    // (A settings-only Session's projection holds its own `ui` type.)
    const requestInputField = !settingsOnly && canonicalProjection?.inputType
      ? (canonicalProjection.mode === 'linear' ? 'lInputType' : 'cInputType') : null;
    /** @param {'cInputType' | 'lInputType'} field */
    const requestInputTypeWins = (field) => field === requestInputField && readSessionInputType(ui[field]).lost;
    /** @param {'cInputType' | 'lInputType'} field @returns {string} */
    const inputTypeOf = (field) => (field === requestInputField && (isUnsetValue(ui[field]) || requestInputTypeWins(field))
      ? /** @type {string} */ (canonicalProjection?.inputType) : readSessionInputType(ui[field]).value);
    const inputTypes = { cInputType: inputTypeOf('cInputType'), lInputType: inputTypeOf('lInputType') };
    // Where the input files are in the file.
    const filesAt = canonicalSession ? 'webFiles.bindings' : 'files';
    const candidateFiles = {
      // `applyFiles` fills the file slots.
      files: /** @type {Record<string, unknown>} */ ({}), linearSeqs: [], circularRecordList: { value: [] },
      circularRecordDiscovery: {}, linearReorderNotice: { value: '' },
      cInputType: { value: inputTypes.cInputType },
      // Read as the drawing reads them (a list of another kind failed the Load).
      linearRecordRows: readLinearRecordRows(linearDraftConfig?.linearRecordLayout?.rows).value,
      linearComparisonPlan: normalizeLinearComparisonPlan(linearDraftConfig?.linearComparisonPlan)
    };
    recordSessionLifecycleEvent('session-candidate-files-start');
    const { collapsedLinearSeqs } = applyFiles(
      canonicalSession ? projectionResult.restoredFiles : data.files,
      {
        adoptCanonicalPayloads: currentSchemaSession, resolveRecordInputs: !settingsOnly,
        targetState: candidateFiles, targetDrawing: candidateFiles, repairs, at: filesAt
      }
    );
    recordSessionLifecycleEvent('session-candidate-files-end');
    const importedResults = canonicalSession
      ? projectionResult.artifactState.results
      : data.results;
    const logicalImportedResults = normalizeLogicalResults(
      (Array.isArray(importedResults) ? importedResults : []).map((result, index) => ({
        name: result?.name || `Result ${index + 1}`,
        content: result?.content || ''
      }))
    );
    const storedEditorState = canonicalSession
      ? projectionResult.artifactState.editorState
      : data.editorState;
    /** @type {ReturnType<typeof featureStateFromCatalog> | null} */
    let currentCatalogFeatureState = null;
    let validatedSessionCatalog = currentSchemaSession
      ? projectionResult?.validatedFeatureCatalog || null
      : null;
    let restoredEditorState = storedEditorState;
    if (currentSchemaSession && validatedSessionCatalog) {
      currentCatalogFeatureState = featureStateFromCatalog(
        validatedSessionCatalog,
        { mode: committedMode }
      );
    }

    const artifactFeatureState = canonicalSession
      ? currentSchemaSession
        ? { ...projectionResult.artifactState.features }
        : cloneJsonData(projectionResult.artifactState.features)
      : {};
    if (currentSchemaSession) {
      [
        'extractedFeatures',
        'biologicalFeatures',
        'featureSelectorSafetyScope',
        'featureRecordIds'
      ].forEach((field) => delete artifactFeatureState[field]);
    }
    let features = splitLegacyFeatureVisibilityRules(canonicalSession
      ? featureEditsOverProjection({ ...(currentCatalogFeatureState || {}), ...artifactFeatureState },
        projectionResult.renderState)
      : (data.features || {}));
    const catalogSequenceSources = currentSchemaSession
      ? (currentCatalogFeatureState?.sequenceSources || [])
      : [];
    recordSessionLifecycleEvent('session-candidate-sequences-start');
    // E1: each Result set reads its sources in its own mode.
    const setSequenceSourceOptions = {
      olderSession: sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION,
      settingsOnly,
      cInputType: candidateFiles.cInputType.value,
      lInputType: inputTypes.lInputType,
      files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs,
      circularConservation: circularDraftConfig?.circularConservation || {}
    };
    const { restored: restoredFileSequenceSources, error: currentRecoveryError } = await restoreLoadedSetSequenceSources({
      ...setSequenceSourceOptions,
      mode: generatedMode,
      catalog: currentSchemaSession ? validatedSessionCatalog : null,
      renderRequest: data.renderRequest
    });
    // OV-712: listed, and reported as `degradedRecovery`.
    if (currentRecoveryError) {
      repairs.push(loadRepair({ code: 'match-sequences-unrestored', path: currentSchemaSession ? 'resources' : 'files', effect: 'dropped' }));
    }

    recordSessionLifecycleEvent('session-candidate-sequences-end');
    // A saved Result is laid out at its committed legend and title sides.
    const committedLayout = canonicalProjection?.layoutPreferences
      ? resolveActiveLayoutPreference(
          canonicalProjection.layoutPreferences,
          canonicalProjection.mode,
          Boolean(canonicalProjection.config?.form?.multi_record_canvas)
        )
      : null;
    // The composition and stroke owners transform an older Session's Result
    // through the composition root's port (R13); this module passes data only.
    /** @type {LegacyResultSvgTransform} */
    const transformLegacyResultSvg = options.transformLegacyResultSvg;
    // A Result of `mode`; `composition` is the legacy composition of the
    // committed set (null for the other set, which has none).
    // A Result's set is the top-level one (`results`) or the other mode's.
    // A Result Load cannot admit or compose drops its set's Results (A20,
    // OV-789): the first reading names the set, and Load reads the file again.
    /** @template R @param {ResultSetPath} at @param {() => R} admit @returns {R} */
    const admittedOrUncomposed = (at, admit) => {
      try {
        return admit();
      } catch (error) {
        if (error instanceof UncomposedResultSet || uncomposedSets.includes(at)) throw error;
        throw new UncomposedResultSet(at, error);
      }
    };
    /** @param {'circular' | 'linear'} mode @param {LegacyResultSvgData['composition']} composition */
    const loadedResultTransform = (mode, composition) => (/** @type {SVGSVGElement} */ svg) => {
      /** @type {ResultSetPath} */
      const at = mode === generatedMode ? '' : 'otherModeResult.';
      return admittedOrUncomposed(at, () => {
        const legendGroupsChanged = composition ? normalizeLegacyLegendEntryGroups(svg) : false;
        const legacyResultChanged = transformLegacyResultSvg(svg, { composition, mode, repairs, path: `${at}results` });
        return legendGroupsChanged || legacyResultChanged;
      });
    };
    /** @param {'circular' | 'linear'} mode @param {LegacyResultSvgData['composition']} composition @returns {LegacyResultNormalization | null} */
    const loadedResultNormalization = (mode, composition) => (transformLegacyResultSvg
      ? {
          transform: loadedResultTransform(mode, composition),
          applies: (content) => transformLegacyResultSvg.appliesToContent(content, mode, Boolean(composition))
        }
      : null);
    // The side a Session 31+ stored for the mode, or that a Session 27-30
    // stored in its form, read through the Result's mode (OV-671); a side
    // outside it is listed when a Result is drawn with it (BA-D26). The form
    // of a Session 27-30 saved after a mode switch is the shown mode's: its
    // Result was drawn with `ui.generatedLegendPosition` (OV-677), read as a
    // drawn side (blank is none) and listed once, above, and its mode's title
    // side, which 0.13.0 saved in `ui` (OV-969).
    const legacyDrawnSide = formIsResultMode
      ? readDrawnLegendSideOfMode(restoredConfig?.form?.legend, generatedMode, 'config.form.legend')
      : { value: drawnLegendSideOfMode(ui.generatedLegendPosition, generatedMode), repairs: [] };
    /** @type {string} */
    const savedTitleSideOfResultMode = generatedMode === 'circular'
      ? readGeneratedCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition, ui).value
      : ui.linearPlotTitlePosition ? normalizeLinearPlotTitlePosition(ui.linearPlotTitlePosition) : 'none';
    if (!committedLayout && logicalImportedResults.length > 0) repairs.push(...legacyDrawnSide.repairs);
    /** @type {LegacyResultSvgData['composition']} */
    const committedLegacyComposition = {
      legendSide: committedLayout?.legend ?? legacyDrawnSide.value,
      titleSide: committedLayout?.plotTitlePosition
        || (formIsResultMode ? restoredConfig?.adv?.plot_title_position : savedTitleSideOfResultMode) || 'none',
      userDeltas: {
        primary: ui.diagramOffset ? [ui.diagramOffset.x, ui.diagramOffset.y] : null,
        legend: ui.legendCurrentOffset
          ? [ui.legendCurrentOffset.x, ui.legendCurrentOffset.y]
          : null,
        lengthBar: ui.lengthBarUserOffset
          ? [ui.lengthBarUserOffset.x, ui.lengthBarUserOffset.y]
          : null,
        title: ui.plotTitleUserOffset
          ? [ui.plotTitleUserOffset.x, ui.plotTitleUserOffset.y]
          : null
      }
    };
    const transformRestoredSessionSvg = loadedResultTransform(generatedMode, committedLegacyComposition);

    recordSessionLifecycleEvent('svg-admission-start');
    // The loaded Results carry the feature types of the request that becomes
    // the committed request below (R13). An older Session adopts no catalog,
    // so its Features list lists what each Result renders and reads none.
    const selectedFeatureTypes = currentSchemaSession
      ? adoptedCanonicalSession?.renderRequest?.diagramOptions?.selectedFeaturesSet
      : null;
    const committedSavedEdits = savedResultEdits(
      data, generatedMode, data.editorState, sourceSessionVersion, draftConfigOf(generatedMode)?.rules
    );
    const committedImportedResults = currentSchemaSession && validatedSessionCatalog
      ? admitLoadedSetResults(logicalImportedResults, {
          featureCatalog: validatedSessionCatalog, mode: committedMode, selectedFeatureTypes,
          savedEdits: committedSavedEdits,
          legacy: loadedResultNormalization(generatedMode, committedLegacyComposition)
        })
      : admittedOrUncomposed('', () => admitLegacyImportedResults(
          createLegacyImportResultSource(logicalImportedResults),
          { transformSvg: transformRestoredSessionSvg, selectedFeatureTypes }
        ));
    recordSessionLifecycleEvent('svg-admission-end');
    // E1: a Session with a Result set of each mode builds both sets' artifact
    // slots with one function. It opens on its shown mode when that mode has a
    // Result, otherwise on the mode that has one; the opening slot is installed
    // last, before the preview mounts, and the other slot waits in its mode.
    const otherModeCatalogAdmitted = currentSchemaSession ? otherModeCatalog : null;
    const otherModeResult = otherModeCatalogAdmitted ? data.otherModeResult : null;
    /** @type {'circular' | 'linear'} */
    const otherSetMode = otherModeResult?.renderRequest?.mode === 'linear' ? 'linear' : 'circular';
    // A projection's mode is checked (`projectCanonicalSessionRequest`) or read (K19).
    /** @type {Set<DiagramMode>} */
    const resultModes = new Set([
      ...(committedImportedResults.length > 0 ? [generatedMode] : []),
      ...(otherModeResult ? [otherSetMode] : [])
    ]);
    const displayMode = resultModes.size === 0 || resultModes.has(candidateMode)
      ? candidateMode
      : [...resultModes][0];
    // Each set's History byte estimate: its share of the decompressed file,
    // split by the size of its Results (a serialization of the other set only
    // to size it costs about 20 ms for a Vibrio-size set).
    /** @param {unknown} setResults */
    const resultCharacters = (setResults) => (Array.isArray(setResults) ? setResults : [])
      .reduce((sum, result) => sum + String(result?.content || '').length, 0);
    const otherResultCharacters = otherModeResult ? resultCharacters(otherModeResult.results) : 0;
    const otherSetCharacters = otherModeResult ? Math.round(candidate.characters * otherResultCharacters
      / Math.max(1, otherResultCharacters + resultCharacters(data.results))) : 0;
    const topSetCharacters = Math.max(0, candidate.characters - otherSetCharacters);
    const slotOptions = {
      resources: data.resources, resourceTable: currentResourceTable, webFiles: data.webFiles
    };
    const otherSavedEdits = otherModeResult
      ? savedResultEdits(data, otherSetMode, otherModeResult.editorState, sourceSessionVersion, draftConfigOf(otherSetMode)?.rules)
      : null;
    // The rule rows of each drawing whose keys the reading of its saved
    // captions changes (F2): its Result set's, else over no generated rows.
    /** @param {'circular' | 'linear'} mode @returns {ReadonlyMap<PythonLegendKey, PythonLegendKey>} */
    const legendRekeysOf = (mode) => ((mode === generatedMode ? committedSavedEdits : mode === otherSetMode ? otherSavedEdits : null)
      ?.legendRekeys ?? savedLegendRekeys(draftConfigOf(mode)?.rules, { originalOrder: [], originalColors: {} }));
    const otherSequenceSources = otherModeResult && otherModeCatalogAdmitted
      ? await restoreLoadedSetSequenceSources({
          ...setSequenceSourceOptions,
          mode: otherSetMode,
          catalog: otherModeCatalogAdmitted,
          renderRequest: otherModeResult.renderRequest
        })
      : null;
    if (otherSequenceSources?.error) {
      repairs.push(loadRepair({ code: 'match-sequences-unrestored', path: 'resources', effect: 'dropped' }));
    }
    const otherArtifactSlot = otherModeResult && otherModeCatalogAdmitted && otherSequenceSources
      ? buildLoadedArtifactSlot(otherModeResult, {
          ...slotOptions,
          featureCatalog: otherModeCatalogAdmitted,
          savedEdits: otherSavedEdits,
          legacy: loadedResultNormalization(otherSetMode, null),
          repairs,
          restoredSequenceSources: otherSequenceSources.restored,
          retainedBytes: otherSetCharacters * 2
        })
      : null;
    const opensOtherSet = Boolean(otherArtifactSlot) && displayMode === otherSetMode;
    // The Linear set's committed request owns the Linear draft's alignment plan
    // and record translations, wherever the set sits; both are admitted as the
    // request projection admits the top-level set's.
    const linearSetLayout = otherModeResult && otherSetMode === 'linear'
      ? canonicalLinearRecordLayout(otherModeResult.renderRequest) : null;

    /** @type {Record<string, any> | null} */
    let legacyFeatureRecoveryPlan = null;
    const legacyFeatureSnapshot = {
      mode: generatedMode, cInputType: candidateFiles.cInputType.value,
      lInputType: inputTypes.lInputType, files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs, results: committedImportedResults,
      selectedResultIndex: ui.selectedResultIndex || 0,
      featureState: features, editorState: restoredEditorState
    };
    if (sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
      try {
        legacyFeatureRecoveryPlan = await buildSessionFeatureRecoveryPlan(/** @type {any} */ ({
          snapshot: legacyFeatureSnapshot,
          // After `splitLegacyFeatureVisibilityRules`, a rule list is the manual
          // rules or none (a value of another kind serialized as no rule).
          featureVisibilityTsv: serializeFeatureVisibilityRules(
            Array.isArray(features.featureVisibilityManualRules) ? features.featureVisibilityManualRules : []
          )
        }));
      } catch (error) {
        legacyFeatureRecoveryPlan = { status: 'failed', warning: normalizeCaughtError(error).summary };
      }
      // The Features card says why after the Load (`applySessionFeatureRecoveryPlan`).
      if (legacyFeatureRecoveryPlan?.status === 'unrecoverable' || legacyFeatureRecoveryPlan?.status === 'failed') {
        repairs.push(loadRepair({ code: 'feature-recovery-failed', path: 'features', effect: 'dropped' }));
      }
    }
    // Session 46 keys per-feature edits by source identity: an older Session's
    // saved rendered-ID edits are mapped through its saved catalog, or without
    // one through its sources read again with its crops and orientations (the
    // saved metadata when they cannot be read), once (design Q4 4.3).
    let migratedAnnotationTargetCount = 0;
    if (sourceSessionVersion < MODE_SCOPED_SESSION_VERSION) {
      const recovered = legacyFeatureRecoveryPlan?.recoveredFeatureState;
      const sessionMode = data.renderRequest?.mode || generatedMode;
      const hashColorRules = sourceSessionVersion <= DRAWN_HASH_RULE_SESSION_MAX_VERSION
        && hasHashColorRules(restoredConfig?.rules);
      const sourceFeatures = !validatedSessionCatalog && (hasRenderedIdFeatureEdits(features) || hashColorRules)
        ? await extractSessionSourceFeatures(/** @type {any} */ ({ snapshot: legacyFeatureSnapshot }))
        : null;
      const legacy = validatedSessionCatalog ? null : {
        records: data.renderRequest?.records,
        features: [
          sourceFeatures?.extractedFeatures, recovered?.biologicalFeatures, recovered?.extractedFeatures,
          features.biologicalFeatures, features.extractedFeatures
        ].find((list) => Array.isArray(list) && list.length > 0) || [],
        biologicalFeatures: sourceFeatures?.biologicalFeatures || []
      };
      const migration = migrateSessionFeatureEdits({
        features, mode: sessionMode, catalog: validatedSessionCatalog, legacy
      });
      features = migration.features;
      if (migration.droppedCount > 0) {
        repairs.push(loadRepair({ code: 'feature-edit-unmatched', path: 'features', effect: 'dropped', count: migration.droppedCount }));
      }
      if (migration.narrowedVisibilityCount > 0) {
        repairs.push(loadRepair({
          code: 'feature-visibility-narrowed', path: 'features.featureVisibilityOverrides', effect: 'changed',
          count: migration.narrowedVisibilityCount
        }));
      }
      // S6: a `hash` color rule of a Session 44 or older names its features by
      // their source hash (OV-401, release D-39). A Session before 31 has no
      // request: its Linear cards give the records.
      if (hashColorRules && restoredConfig) {
        const ruleRecords = data.renderRequest?.records
          ?? (generatedMode === 'linear' ? legacyLinearRequestRecords(candidateFiles.linearSeqs) : undefined);
        const hashRules = migrateSessionHashRules({
          rules: restoredConfig.rules,
          mode: sessionMode,
          catalog: validatedSessionCatalog,
          legacy: legacy && { ...legacy, records: ruleRecords },
          records: ruleRecords
        });
        restoredConfig.rules = hashRules.rules;
        if (hashRules.unmappedCount > 0) {
          repairs.push(loadRepair({ code: 'hash-selector-unmapped', path: 'config.rules', effect: 'changed', count: hashRules.unmappedCount }));
        }
      }
      // R-7: an annotation's `hash=` target moves to its source feature only
      // where the saved catalog makes the figure certain.
      if (restoredConfig) {
        const annotationMigration = migrateSessionAnnotationTargets({
          annotationSets: restoredConfig.annotationSets,
          mode: sessionMode,
          catalog: validatedSessionCatalog,
          records: data.renderRequest?.records
        });
        migratedAnnotationTargetCount = annotationMigration.migratedCount;
        if (migratedAnnotationTargetCount > 0) restoredConfig.annotationSets = annotationMigration.annotationSets;
      }
      if (recovered) {
        const recoveredWithoutRenderedIdEdits = { ...recovered };
        RENDERED_ID_FEATURE_EDIT_FIELDS.forEach((field) => delete recoveredWithoutRenderedIdEdits[field]);
        // `recovered` came from legacyFeatureRecoveryPlan.recoveredFeatureState, so the plan is set.
        // Its edit rows are the committed mode's drawing's (R1-3).
        /** @type {Record<string, any>} */ (legacyFeatureRecoveryPlan).recoveredFeatureState = {
          ...recoveredWithoutRenderedIdEdits,
          featureOverrides: unscopedDraftRows(features.featureOverrides, generatedMode),
          labelOverrideRows: features.labelOverrideRows
        };
      }
    }
    // Each drawing's draft (PD-OI-086, plan 4.1). A Session 46 slice lies over
    // its mode's committed projection. The flat draft of a Session 27-44, as
    // Load has read it so far, is split by the registry (plan 4.2); a Session
    // written without a draft (CLI, Python) splits only what its migrations
    // wrote, and its committed mode keeps its request's projection.
    const hasFlatDraft = !modeScopedSession && isPlainObject(data.config);
    /** @type {Record<'circular' | 'linear', Record<string, any>>} */
    let modeSlices = { circular: {}, linear: {} };
    /** @type {Record<string, any>} */
    let appUi = modeScopedSession && isPlainObject(data.ui) ? data.ui : {};
    let sessionCliOptions = modeScopedSession ? data.cliOptions : undefined;
    // The stored draft a Session 27-44's drawings both read.
    /** @type {unknown} */
    let flatDraft = null;
    if (modeScopedSession) {
      const modes = /** @type {UnvalidatedFields<'circular' | 'linear'>} */ (data.modes);
      modeSlices = {
        circular: isPlainObject(modes?.circular) ? modes.circular : {},
        linear: isPlainObject(modes?.linear) ? modes.linear : {}
      };
    } else {
      const flatConfig = hasFlatDraft
        ? restoredConfig
        : migratedAnnotationTargetCount > 0 && restoredConfig ? { annotationSets: restoredConfig.annotationSets } : {};
      flatDraft = hasFlatDraft ? restoredConfig : null;
      const legend = isPlainObject(restoredEditorState?.legend) ? restoredEditorState.legend : {};
      // The split writes what the Session saved and its migrations wrote
      // (plan 4.2), not the empty defaults Load reads in their place.
      const savedDraft = (/** @type {unknown} */ value) => isPlainObject(value);
      const storedEditorEdits = /** @type {UnvalidatedFields<'legend' | 'featureStrokes'>} */ (data.editorState);
      /**
       * @param {Record<string, any> | null | undefined} value
       * @param {unknown} saved
       * @param {string[]} keys
       */
      const withoutUnsavedEmpty = (value, saved, keys) => {
        if (!value || !isPlainObject(value)) return value;
        const draft = /** @type {Record<string, any>} */ (value);
        const unsaved = keys.filter((key) => {
          if (!Object.hasOwn(draft, key) || (isPlainObject(saved) && Object.hasOwn(/** @type {Record<string, any>} */ (saved), key))) return false;
          const field = draft[key];
          return isPlainObject(field) ? Object.keys(field).length === 0
            : Array.isArray(field) ? field.length === 0 : field === 0;
        });
        return Object.fromEntries(Object.entries(draft).filter(([key]) => !unsaved.includes(key)));
      };
      const split = splitDraftIntoModes({
        // A Session 27-44 draft never saved `fileImports` (Session 48).
        config: withoutUnsavedEmpty(flatConfig, data.config, ['unmanagedConfigOverrides', 'fileImports']) || {},
        // With a catalog, the feature-edit migration writes `featureOverrides`.
        ...(savedDraft(data.features) ? { features: withoutUnsavedEmpty(features, data.features, [
          'featureColorOverrides', 'featureVisibilityManualRules', 'labelOverrideRows', 'labelTextBulkOverrides',
          'selectedFeatureRecordIdx', ...(validatedSessionCatalog ? [] : ['featureOverrides'])
        ]) } : {}),
        editorState: {
          ...(savedDraft(storedEditorEdits?.legend) ? { legend } : {}),
          ...(savedDraft(storedEditorEdits?.featureStrokes)
            ? { featureStrokes: restoredEditorState?.featureStrokes || { overrides: {} } } : {})
        },
        // A saved layout goes to its slots (registry `layout`); legacy layout
        // fields are restored below, as before.
        ui,
        renderRequest: data.renderRequest
      }, {
        committedMode: generatedMode,
        depthSources: sessionDepthSourceWidths({
          c_depth: candidateFiles.files.c_depth, linearSeqs: candidateFiles.linearSeqs
        }),
        paletteColors: paletteColorsFromDefinitions(String(flatConfig?.palette || '').trim() || 'default')
      });
      // The shared migration vectors read the split as Load made it
      // (tests/web/mode-split-vectors.test.mjs).
      recordSessionLifecycleEvent('session-draft-split', { split });
      modeSlices = split.modes || modeSlices;
      appUi = split.ui || {};
      sessionCliOptions = split.cliOptions;
    }
    // The draft each drawing installs: its slice over its committed projection.
    /** @type {Record<'circular' | 'linear', Record<string, any>>} */
    const modeConfigs = {
      circular: {}, linear: {}
    };
    // Each drawing's editor edits, read here so what they lose is listed.
    /** @type {Record<'circular' | 'linear', ReturnType<typeof normalizeEditorStateData>>} */
    const modeEditorStates = {
      circular: normalizeEditorStateData({}), linear: normalizeEditorStateData({})
    };
    // The layout of the committed request (a Session 46 slice without its
    // own slot takes it).
    const projectedLayout = canonicalSession ? canonicalProjection?.layoutPreferences : null;
    /** @type {Record<'circular' | 'linear', SettingsContainers>} */
    const installedSettings = { circular: defaultSettingsContainers('circular'), linear: defaultSettingsContainers('linear') };
    // The legend positions of the drafts' layouts: a Session 46 slice holds
    // its mode's slot; a Session 27-44 keeps one set in `ui` (OV-652).
    if (modeScopedSession) {
      SLICE_MODES.forEach((mode) => repairs.push(...storedLegendPositionRepairs(modeSlices[mode].ui, `modes.${mode}.ui`, {
        slotOf: mode, projected: mode === generatedMode && Boolean(projectedLayout)
      })));
    } else {
      repairs.push(...storedLegendPositionRepairs(ui, 'ui'));
    }
    for (const mode of SLICE_MODES) {
      const storedEditorState = modeSlices[mode].editorState;
      const editorPath = modeScopedSession ? `modes.${mode}.editorState` : 'editorState';
      const editorState = normalizeEditorStateData(isPlainObject(storedEditorState) ? storedEditorState : {}, {
        repairs, path: editorPath, sessionVersion: sourceSessionVersion
      });
      // The rule rows whose key the reading of the saved captions changes
      // take the Legend records stored under the old key (F2); an edit the
      // new key already has is kept, and the moved one is listed (B39).
      /** @type {DroppedLegendRecord[]} */
      const droppedLegendRecords = [];
      modeEditorStates[mode] = {
        ...editorState, legend: rekeyLegendRecords(editorState.legend, legendRekeysOf(mode), droppedLegendRecords)
      };
      droppedLegendRecords.forEach(({ collection, from, to }) => repairs.push(loadRepair({
        code: 'legend-edit-rekey-dropped', path: `${editorPath}.legend.${collection}`, effect: 'dropped', label: from, used: to
      })));
      const committed = mode === generatedMode;
      modeConfigs[mode] = modeScopedSession
        ? (committed && committedSliceProjected ? restoredConfig : overlayModeSliceConfig(null, modeSlices[mode].config))
        : overlayModeSliceConfig(committed && !hasFlatDraft ? restoredConfig : null, modeSlices[mode].config);
      // A caption an older build saved with its own rule's color ID is read without it (PD-OI-042 scenario revision 2).
      if (Array.isArray(modeConfigs[mode].rules)) modeConfigs[mode].rules = savedSpecificRules(modeConfigs[mode].rules);
      if (mode === 'linear') modeConfigs.linear = withBoundLinearRecordsOff(modeConfigs.linear, candidateFiles.linearSeqs);
      if (mode === 'circular') {
        modeConfigs.circular = withoutUnreadConservationSeries(modeConfigs.circular, candidateFiles.files, {
          seriesAt: `${modeScopedSession ? 'modes.circular.config' : 'config'}.circularConservation.series`,
          filesAt,
          stored: canonicalSession ? null : /** @type {UnvalidatedFields<'c_conservation_blasts'>} */ (data.files)?.c_conservation_blasts
        }, repairs);
      }
      // The form's legend position is the drawing's layout slot (state.js), which reads it in its mode.
      if (isPlainObject(modeConfigs[mode].form)) {
        repairs.push(...readLegendPositionOfMode(modeConfigs[mode].form.legend, mode,
          `${modeScopedSession ? `modes.${mode}.config` : 'config'}.form.legend`).repairs);
      }
      // What the settings' readers cannot restore of the stored draft (B37,
      // B43), installed as `applyModeSliceData` installs it over the drawing
      // Load resets; and of the drawing's `ui` and feature edits.
      installedSettings[mode] = defaultSettingsContainers(mode);
      // What `applyModeSliceData` installs.
      const installedConfig = overlayModeSliceConfig(null, modeConfigs[mode]);
      repairs.push(...configFieldRepairs(installedSettings[mode], installedConfig,
        modeScopedSession ? `modes.${mode}.config` : 'config', {
          stored: modeSlices[mode].config, mode, resolveTrackPlacements: !settingsOnly && !isStoredSlice(mode)
        }));
      repairs.push(...unreadFieldRepairs(modeScopedSession ? modeSlices[mode].config : flatDraft,
        modeScopedSession ? `modes.${mode}.config` : 'config', mode, { slice: modeScopedSession }));
      // (A Session 27-44 keeps them in its `ui` and `features`.)
      repairs.push(...sliceFieldRepairs(modeSlices[mode], modeScopedSession ? `modes.${mode}` : '', mode));
      // A saved unmanaged override is checked against its own drawing's mode (PD-OI-009, OV-106).
      const savedOverrides = modeSlices[mode].config?.unmanagedConfigOverrides;
      if (isPlainObject(savedOverrides) && Object.keys(savedOverrides).length > 0) {
        modeConfigs[mode].unmanagedConfigOverrides = await readStoredUnmanagedConfigOverrides({
          mode, configOverrides: savedOverrides, source: 'draft'
        }, `${modeScopedSession ? `modes.${mode}.config` : 'config'}.unmanagedConfigOverrides`, repairs);
      }
    }
    // What each drawing's Depth fit drops (B45): it runs on the restored files.
    if (!settingsOnly) {
      const linearRows = candidateFiles.linearSeqs.map((/** @type {LinearSeq} */ seq) => depthFileSlotsFromValue(seq.depth));
      // `applyFiles` restored the bound Circular Depth slots (by record, or one row).
      const circularFiles = /** @type {DepthFileSlots | DepthFileSlots[]} */ (/** @type {{ c_depth?: unknown }} */ (candidateFiles.files).c_depth);
      SLICE_MODES.forEach((mode) => repairs.push(...depthFitRepairs(mode, installedSettings[mode], modeSlices[mode].config?.adv, {
        at: `${modeScopedSession ? `modes.${mode}.config` : 'config'}.adv`,
        resolveTrackPlacements: !isStoredSlice(mode),
        circularFiles,
        linearRows
      })));
    }
    // The app settings Load reads from `ui` (B36, B43), and a drawing's
    // pending palette, which Instant Preview drops (B44): the drafts apply in
    // SLICE_MODES order, each over the Instant Preview the one before left.
    // A boolean Instant Preview of an older draft wins over `ui`'s.
    const draftInstantPreview = (appUi !== ui && typeof appUi.paletteInstantPreviewEnabled === 'boolean')
      || SLICE_MODES.some((mode) => typeof modeConfigs[mode].paletteInstantPreviewEnabled === 'boolean');
    repairs.push(
      ...sessionUiFieldRepairs({
        ...ui,
        ...(draftInstantPreview ? { paletteInstantPreviewEnabled: undefined } : {}),
        // The request's input type wins over a stored one of no type.
        ...Object.fromEntries(/** @type {const} */ (['cInputType', 'lInputType']).filter(requestInputTypeWins).map((field) => [field, undefined]))
      }, 'ui'),
      ...losatExecutionRepairs(appUi.losatExecution, modeScopedSession ? 'ui.losatExecution' : 'config.losat', { own: modeScopedSession }),
      ...valueRepairs(readRichFeaturePopup, appUi.richFeaturePopup, modeScopedSession ? 'ui.richFeaturePopup' : 'config.adv.rich_feature_popup')
    );
    let instantPreview = Boolean(typeof appUi.paletteInstantPreviewEnabled === 'boolean'
      ? appUi.paletteInstantPreviewEnabled : ui.paletteInstantPreviewEnabled);
    for (const mode of SLICE_MODES) {
      if (typeof modeConfigs[mode].paletteInstantPreviewEnabled === 'boolean') instantPreview = modeConfigs[mode].paletteInstantPreviewEnabled;
      const pendingPalette = pendingPaletteName(modeSlices[mode].ui);
      const pendingPath = `${modeScopedSession ? `modes.${mode}.ui` : 'ui'}`;
      if (pendingPalette && instantPreview) {
        repairs.push(loadRepair({ code: 'pending-palette-dropped', path: `${pendingPath}.pendingPaletteName`, effect: 'dropped', label: pendingPalette }));
      } else if (pendingPalette) {
        // The pending palette's colors apply with it.
        repairs.push(...valueRepairs(readPaletteColors, modeSlices[mode].ui.pendingPaletteColors, `${pendingPath}.pendingPaletteColors`));
      }
    }
    // The similarity groups' state Load installs; a name or description edit
    // of a group it does not hold is dropped and listed (B47).
    const loadedOrthogroupState = canonicalSession
      ? {
          ...projectionResult.artifactState.orthogroupState,
          ...(modeScopedSession ? {
            orthogroupNameOverrides: modeConfigs.linear.webEdits?.orthogroupNameOverrides || {},
            orthogroupDescriptionOverrides: modeConfigs.linear.webEdits?.orthogroupDescriptionOverrides || {}
          } : {}),
          ...(currentCatalogFeatureState
            ? { groups: currentCatalogFeatureState.orthogroups }
            : {})
        }
      : data.orthogroupState && typeof data.orthogroupState === 'object'
        ? data.orthogroupState
      : {
          groups: [],
          selectedOrthogroupId: features.selectedOrthogroupId,
          selectedOrthogroupAlignmentFeature: features.selectedOrthogroupAlignmentFeature,
          orthogroupNameOverrides:
            features.orthogroupNameOverrides ||
            /** @type {UnvalidatedFields<'webEdits', UnvalidatedFields<'orthogroupNameOverrides'>>} */ (data.config)
              ?.webEdits?.orthogroupNameOverrides ||
            {},
          orthogroupDescriptionOverrides:
            features.orthogroupDescriptionOverrides ||
            /** @type {UnvalidatedFields<'webEdits', UnvalidatedFields<'orthogroupDescriptionOverrides'>>} */ (data.config)
              ?.webEdits?.orthogroupDescriptionOverrides ||
            {}
      };
    {
      const groupIds = new Set(orthogroupIdsOf(loadedOrthogroupState.groups));
      const editsPath = canonicalSession
        ? (modeScopedSession ? 'modes.linear.config.webEdits' : 'orthogroupState')
        : (isPlainObject(data.orthogroupState) ? 'orthogroupState' : 'features');
      /** @type {['orthogroupNameOverrides' | 'orthogroupDescriptionOverrides', unknown][]} */ ([
        ['orthogroupNameOverrides', loadedOrthogroupState.orthogroupNameOverrides],
        ['orthogroupDescriptionOverrides', loadedOrthogroupState.orthogroupDescriptionOverrides]
      ]).forEach(([field, edits]) => {
        const { dropped } = orthogroupEditsOf(edits, groupIds);
        if (dropped) repairs.push(loadRepair({ code: 'orthogroup-edit-dropped', path: `${editsPath}.${field}`, effect: 'dropped', count: dropped }));
      });
    }
    // The LOSAT caches and the protein identity manifest Load drops are
    // computed again by the next comparison: no dialog row (Q1, B46).
    const cacheSource = canonicalSession ? projectionResult.artifactState : data;
    const loadedProteinManifest = cacheSource.proteinIdentityManifest;
    const loadedProteinManifestValid = validateProteinIdentityManifest(loadedProteinManifest);
    /** @type {[string, number][]} */
    const droppedCacheEntries = [
      ['losatCache.entries', (Array.isArray(cacheSource.losatCache?.entries) ? cacheSource.losatCache.entries : [])
        .filter((/** @type {unknown} */ entry) => savedLosatCacheEntryKind(entry) === 'dropped').length],
      ['losatDerivedCache.entries', readLosatDerivedCache(cacheSource.losatDerivedCache?.entries, { adoptCurrent: true }).dropped],
      // A manifest the file holds (the projection reads an absent one as empty).
      ['proteinIdentityManifest', !loadedProteinManifestValid && data.proteinIdentityManifest != null
        && !(isPlainObject(data.proteinIdentityManifest) && !Object.keys(data.proteinIdentityManifest).length) ? 1 : 0]
    ];
    droppedCacheEntries.forEach(([path, count]) => {
      if (count) repairs.push(loadRepair({ code: 'cache-entry-dropped', path, effect: 'dropped', count }));
    });
    // The top-level set, built like the other set, waits in its mode when the
    // Session opens on the other set's mode. The stash starts again with the
    // loaded Session. Its editor state is the one read above.
    const topArtifactSlot = opensOtherSet && validatedSessionCatalog
      ? buildLoadedArtifactSlot({
          renderRequest: data.renderRequest,
          results: data.results,
          editorState: data.editorState,
          ui,
          runMetadata: projectionResult?.artifactState?.runMetadata
        }, {
          ...slotOptions,
          featureCatalog: validatedSessionCatalog,
          admittedResults: committedImportedResults,
          restoredSequenceSources: restoredFileSequenceSources,
          retainedBytes: topSetCharacters * 2,
          repairs: null
        })
      : null;
    // The saved Result's editor edits: a Session 39 or older is read here,
    // so what its reading loses is listed.
    const artifactEditorState = currentSchemaSession ? restoredEditorState : normalizeEditorStateData(restoredEditorState, {
      repairs, path: 'editorState', sessionVersion: sourceSessionVersion
    });
    recordSessionLifecycleEvent('session-candidate-prepared');
    // BA-D25: what the candidate could not restore is listed before anything
    // changes; Cancel leaves the settings, Results, inputs, History and screen
    // as they are.
    Object.freeze(repairs);
    // R4: a file of which Load restores nothing is refused (every producer ran; nothing is committed).
    if (restoresNothing(storedDocument, repairs)) throw SESSION_FIELDS_ERROR();
    // One record per place in the file and reading (OV-765).
    const listedRepairs = Object.freeze(distinctLoadRepairs(repairs));
    const repairDialog = describeLoadRepairs(listedRepairs);
    if (repairDialog) {
      const confirmLoadRepairs = /** @type {ConfirmLoadRepairs} */ (options.confirmLoadRepairs);
      if (!await confirmLoadRepairs(repairDialog) || !options.isCurrent()) return { status: 'canceled' };
    }
    // Custom fonts (receipt 〔読み込むとき〕): a font either drawing names that
    // this page cannot draw asks before anything is replaced; Cancel loads nothing.
    // It asks after the repair list, so a Cancel there loads no font either.
    /** @type {SessionFontConfirmation | undefined} */
    const confirmSessionFonts = options.confirmSessionFonts;
    if (confirmSessionFonts && !await confirmSessionFonts(SLICE_MODES.map((mode) => modeConfigs[mode]?.adv?.font_family))) {
      return { status: 'canceled' };
    }
    // The admitted candidate is still private. Let the browser handle input
    // after SVG sanitation, before the one atomic live-state transaction.
    await new Promise((resolve) => setTimeout(resolve, 0));
    if (!options.isCurrent()) return { status: 'canceled' };
    rollbackSnapshot = captureSessionImportSnapshot();
    if (typeof rollbackStateExtension?.capture === 'function') {
      rollbackExtensionSnapshot = rollbackStateExtension.capture();
    }
    commitStarted = true;
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    state.sessionResourceDiscoveryDeferred.value = currentSchemaSession;
    if (currentSchemaSession) {
      setCommittedCanonicalSession(adoptedCanonicalSession);
      activeSessionResourceTable = currentResourceTable;
    }
    state.sessionTitle.value = canonicalSession
      ? projectionResult.documentMetadata.title
      : (typeof data.title === 'string' ? data.title : '');
    if (ui.mode) state.mode.value = ui.mode;
    state.cInputType.value = inputTypes.cInputType;
    state.lInputType.value = inputTypes.lInputType;
    state.downloadDpi.value = readDownloadDpi(ui.downloadDpi).value;
    // Suppressed mode watchers observe only the complete admitted document.
    state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
    // A Session 27-44 draft's boolean Instant Preview wins over its `ui` (plan 4.2).
    state.paletteInstantPreviewEnabled.value = Boolean(
      typeof appUi.paletteInstantPreviewEnabled === 'boolean'
        ? appUi.paletteInstantPreviewEnabled : ui.paletteInstantPreviewEnabled
    );
    state.featurePanelTab.value = readFeaturePanelTab(ui.featurePanelTab).value;
    state.generatedMode.value = generatedMode;
    if (generatedLegend) state.generatedLegendPosition.value = generatedLegend.value;
    state.generatedMultiRecordCanvas.value = Boolean(ui.generatedMultiRecordCanvas);
    state.generatedCircularPlotTitlePosition.value = readGeneratedCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition, ui).value;

    state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(
      state.drawings.circular,
      modeConfigs.circular.form
    );
    // Each drawing over its mode's defaults (resetSessionBaseline); the shown
    // drawing's Canvas padding follows once its Result is mounted. A Session 46
    // slice without its mode's request projection is stored settings: its track
    // stacks install as saved, unset axis indexes included (R11).
    for (const mode of SLICE_MODES) {
      // The committed mode's slice takes its request's per-feature edits where
      // it has none (OV-220).
      const features = mode === generatedMode && canonicalSession && !settingsOnly
        ? { features: featureEditsOverProjection(modeSlices[mode].features, projectionResult.renderState) } : {};
      applyModeSliceData(state.drawings[mode], mode, {
        ...modeSlices[mode], ...features, config: modeConfigs[mode], editorState: modeEditorStates[mode]
      }, {
        resolveTrackPlacements: !settingsOnly && !isStoredSlice(mode),
        applyCanvasPadding: mode !== displayMode
      });
    }
    // App-level settings: LOSAT execution (after the drafts, whose request
    // projection may name a thread count), the rich popup, CLI provenance.
    if (isPlainObject(appUi.losatExecution)) applyLosatExecutionData(appUi.losatExecution);
    state.richFeaturePopup.value = readRichFeaturePopup(appUi.richFeaturePopup).value;
    preservedCliOptions = isPlainObject(sessionCliOptions) ? cloneJsonData(sessionCliOptions) : null;
    state.mode.value = displayMode;
    const shownDrawing = state.drawings[displayMode];
    const canonicalLinearLayout = linearSetLayout || canonicalProjection?.config?.linearRecordLayout;
    if (canonicalLinearLayout && state.linearRecordTranslations) {
      state.linearRecordTranslations.value = cloneJsonData(
        canonicalLinearLayout.recordTranslations || []
      );
    }
    if (canonicalLinearLayout && state.similarityAlignmentPlan) {
      state.similarityAlignmentPlan.value = canonicalLinearLayout.similarityAlignment
        ? cloneJsonData(canonicalLinearLayout.similarityAlignment)
        : null;
    }
    restoreAppliedPaletteFromSession(shownDrawing, ui);
    // Layout preferences: a Session 46 slice holds its mode's slot (the
    // committed mode's request sets it when the slice has none); a Session
    // 27-44 keeps one set, read once, which both drawings take.
    if (modeScopedSession) {
      if (projectedLayout && !isPlainObject(modeSlices[generatedMode].ui?.layoutPreferences)) {
        const drawingOfCommit = state.drawings[generatedMode];
        replaceLayoutPreferences(drawingOfCommit.layoutPreferences, {
          ...cloneJsonData(drawingOfCommit.layoutPreferences),
          [generatedMode]: normalizeLayoutPreferences(projectedLayout)[generatedMode]
        });
      }
    } else {
      SLICE_MODES.forEach((mode) => restoreLayoutPreferences(state.drawings[mode], ui, { projected: projectedLayout }));
    }

    restoreLiveFileState(state.drawings.linear, {
      files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs,
      circularRecordList: candidateFiles.circularRecordList.value,
      circularRecordDiscovery: candidateFiles.circularRecordDiscovery,
      linearRecordRows: candidateFiles.linearRecordRows,
      linearComparisonPlan: candidateFiles.linearComparisonPlan
    });
    // The loaded File cards are new to the user: view state, reset with the document.
    foldLoadedLinearSources(
      state.linearFoldedSources, groupLinearSourceRecords(/** @type {LinearSeq[]} */ (state.linearSeqs)),
      state.lInputType.value
    );
    // The comparisons are the Linear drawing's.
    restoreImportedComparisonIntent(
      state.drawings.linear.importedComparisonIntent,
      comparisonClassification,
      linearDraftConfig?.importedComparisonResolution
    );
    if (!settingsOnly) {
      SLICE_MODES.forEach((mode) => reconcileDepthTrackStateAfterSessionFiles(state.drawings[mode], mode));
    }
    if (canonicalSession) {
      applyLosatCache(
        state.drawings.linear,
        projectionResult.artifactState.losatCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinRawCandidates,
        {
          adoptCurrent: currentSchemaSession,
          validatedManifest: projectionResult.artifactState.proteinIdentityManifest
        }
      );
      applyLosatDerivedCache(
        projectionResult.artifactState.losatDerivedCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinDerivedEvidence,
        { adoptCurrent: currentSchemaSession }
      );
      applyProteinIdentityManifest(loadedProteinManifest, {
        adoptCurrent: currentSchemaSession, valid: loadedProteinManifestValid
      });
    } else {
      const legacyArtifacts = /** @type {UnvalidatedFields<'proteinRawCandidates' | 'proteinDerivedEvidence'>} */ (data.legacyArtifacts);
      applyLosatCache(
        state.drawings.linear,
        /** @type {UnvalidatedFields<'entries'>} */ (data.losatCache)?.entries,
        legacyArtifacts?.proteinRawCandidates
      );
      applyLosatDerivedCache(
        /** @type {UnvalidatedFields<'entries'>} */ (data.losatDerivedCache)?.entries,
        legacyArtifacts?.proteinDerivedEvidence
      );
      applyProteinIdentityManifest(loadedProteinManifest, { valid: loadedProteinManifestValid });
    }
    if (collapsedLinearSeqs) {
      state.losatCacheInfo.value = [];
    }

    state.skipCaptureBaseConfig.value = false;

    // The drawings hold their edits; the Result's catalog features are artifacts.
    applyFeatureArtifactData(features);
    state.selectedFeatureRecordIdx.value = Number.isInteger(modeSlices[displayMode].ui?.selectedFeatureRecordIdx)
      ? modeSlices[displayMode].ui.selectedFeatureRecordIdx : 0;
    if (currentSchemaSession && currentCatalogFeatureState) {
      state.collinearGroups.value = currentCatalogFeatureState.collinearGroups;
      synchronizeRestoredFeatureSummaryStatus({ generationId: 'session-load' });
    }
    state.matchSequenceRegistry?.reset?.([
      ...catalogSequenceSources,
      ...restoredFileSequenceSources
    ]);

    // Group names are the Linear drawing's (Session 46 `webEdits`).
    applyOrthogroupStateData(
      state.drawings.linear,
      loadedOrthogroupState,
      {
        legacyRecords: sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION ? data.renderRequest?.records : null,
        catalogFeatureState: currentCatalogFeatureState
      }
    );
    applyEditorArtifactData({ ...artifactEditorState, legend: rekeyLegendRecords(artifactEditorState.legend, legendRekeysOf(generatedMode)) });
    if (legacyFeatureRecoveryPlan) {
      applySessionFeatureRecoveryPlan(state.drawings[generatedMode], legacyFeatureRecoveryPlan, { generationId: 'session-load' });
    }

    let desiredResultIndex = (
      Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
    )
      ? Math.min(ui.selectedResultIndex, Math.max(0, committedImportedResults.length - 1))
      : 0;
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    state.skipCaptureBaseConfig.value = true;
    recordSessionLifecycleEvent('session-candidate-adopted');
    recordSessionLifecycleEvent('preview-mount-start');
    applyResultsData(committedImportedResults, ui);
    state.annotationWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.annotationWarnings || []
    );
    state.featureIdentityNotices.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.featureIdentityNotices || []
    );
    state.comparisonWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.comparisonWarnings || []
    );
    // A Session does not save the label fit report: a loaded Result has no notice.
    state.labelFitReports.value = [];
    state.trackSlotResolvedGeometry.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.trackSlotGeometry ?? null
    );
    // The composition root keeps the loaded artifacts: the opening one is shown, the stashed one waits.
    /** @type {((slots: { opening: Readonly<ArtifactSlot> | null, stashed: Readonly<ArtifactSlot> | null }) => void) | undefined} */
    const installLoadedArtifactSlots = options.installLoadedArtifactSlots;
    {
      installLoadedArtifactSlots?.(opensOtherSet
        ? { opening: otherArtifactSlot, stashed: topArtifactSlot }
        : { opening: null, stashed: otherArtifactSlot });
      if (opensOtherSet) desiredResultIndex = Number(state.selectedResultIndex.value) || 0;
    }
    const previewReadiness = typeof options?.beforePreviewMount === 'function'
      ? options.beforePreviewMount({
          results: state.results.value,
          resultIndex: desiredResultIndex,
          data,
          ui
        })
      : null;
    recordSessionLifecycleEvent('firstCommittedPreview', {
      resultCount: state.results.value.length
    });
    await nextTick();
    recordSessionLifecycleEvent('preview-mount-end');
    if (previewReadiness?.promise) await previewReadiness.promise;
    else if (previewReadiness?.then) await previewReadiness;

    if (typeof options?.afterLoad === 'function') {
      await options.afterLoad({ data, ui });
    }

    const shownPadding = modeSlices[displayMode].ui?.canvasPadding;
    if (shownPadding !== undefined) Object.assign(shownDrawing.canvasPadding, readCanvasPadding(shownPadding).value);
    if (ui.canvasPan) {
      state.canvasPan.x = ui.canvasPan.x || 0;
      state.canvasPan.y = ui.canvasPan.y || 0;
    }
    if (typeof ui.zoom === 'number') {
      state.zoom.value = ui.zoom;
    }

    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    await nextTick();
    if (!currentSchemaSession) {
      // Only a request Load read and projected is committed: a Session before
      // 31 has none (`readSessionRequests`, OV-862), whatever its file holds (OV-819).
      const stored = canonicalSession ? cloneCanonicalSession(data) : null;
      // OV-1181 B-LOAD: after the stored figure, before History's base (`afterImport`).
      const committed = stored && await convertDetachedProteinResultFrames(
        sourceSessionVersion, stored, canonicalProjectionInput
      );
      if (!committed && logicalImportedResults.length > 0) {
        const { adv } = shownDrawing;
        loadedTrackSettings = Object.freeze({
          mode: displayMode,
          dinucleotide: adv.nt,
          slots: adv[`${displayMode}_track_slots_enabled`] ? cloneJsonData(adv[`${displayMode}_track_slots`]) : null
        });
      }
      setCommittedCanonicalSession(committed);
      activeSessionResourceTable = null;
    } else {
      // OV-1145 B-LOAD: after the stored figure, before History's base (`afterImport`).
      const convertedKeptSlot = await convertLoadedProteinResultFrames(
        sourceSessionVersion, opensOtherSet ? topArtifactSlot : otherArtifactSlot
      );
      if (convertedKeptSlot && options.isCurrent()) installLoadedArtifactSlots?.({ opening: null, stashed: convertedKeptSlot });
    }
    recordSessionLifecycleEvent('interactiveReady', {
      status: 'success',
      degradedRecovery: Boolean(currentRecoveryError)
    });
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    await options.afterImport?.({
      status: 'ok',
      decompressedCharacters: opensOtherSet ? otherSetCharacters : topSetCharacters,
      isCurrent: options.isCurrent
    });
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    return {
      status: 'ok',
      // What the load changed without a loss, for the caller's notice; what
      // it could not restore was listed before the load (`repairs`).
      notices: migratedAnnotationTargetCount > 0 ? [ANNOTATION_TARGET_MIGRATION_NOTICE(migratedAnnotationTargetCount)] : [],
      repairs: listedRepairs,
      data,
      decompressedCharacters: candidate.characters,
      degradedRecovery: Boolean(currentRecoveryError),
      comparisonDisposition: state.drawings.linear.importedComparisonIntent.disposition
    };
  } catch (err) {
    if (err?.name === 'AbortError' && !commitStarted) return { status: 'canceled' };
    // Nothing is committed before the dialog: read the file again without the set (A20).
    if (err instanceof UncomposedResultSet && !commitStarted) {
      return importSessionDocument(e, { ...options, uncomposedSets: [...uncomposedSets, err.at] });
    }
    // The inner operation names the cause; a failure without one is the Load's.
    const caught = normalizeCaughtError(err, { stage: 'request-validation' });
    const error = caught.operation === 'unknown' ? normalizeCaughtError(caught, { operation: 'session-load' }) : caught;
    const currentAlert = state.errorLog.value;
    const canNotify = currentAlert === previousAlert || currentAlert === null;
    if (commitStarted && rollbackSnapshot) {
      try {
        await restoreSessionImportSnapshot(rollbackSnapshot);
        if (typeof rollbackStateExtension?.restore === 'function') {
          await rollbackStateExtension.restore(rollbackExtensionSnapshot);
        }
      } catch (rollbackError) {
        state.generationFailureRecovery.value = 'restore-failed';
      }
    }
    const restoredPreviousAlert = JSON.stringify(normalizeUserFacingError(state.errorLog.value))
      === JSON.stringify(normalizeUserFacingError(previousAlert));
    if (!canNotify) {
      if (restoredPreviousAlert) state.errorLog.value = currentAlert;
      return { status: 'stale' };
    }
    if (state.errorLog.value !== previousAlert && state.errorLog.value !== null && !restoredPreviousAlert) return { status: 'stale' };
    state.errorLog.value = error;
    recordSessionLifecycleEvent('interactiveReady', {
      status: 'error',
      error: error.code
    });
    return { status: 'error', error };
  } finally {
    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    e.target.value = '';
  }
};


/** @type {{ canceled: boolean, promise: Promise<any> | null } | null} */
let sessionSaveInFlight = null;
/** @type {{ canceled: boolean, controller: AbortController } | null} */
let activeSessionImport = null;

export const disposeSessionOperations = () => {
  if (sessionSaveInFlight) sessionSaveInFlight.canceled = true;
  if (activeSessionImport) {
    activeSessionImport.canceled = true;
    activeSessionImport.controller.abort();
  }
  sessionSaveInFlight = null;
  activeSessionImport = null;
  state.sessionSavePending.value = false;
  state.sessionImportPending.value = false;
};

// `options.availability` is the Save and Load availability the composition
// root composes from `sessionOperationAvailability` and the edits still
// applying (R13); `options.recordDisplayRows` are the record display rows of
// both drawings' sources, on which Save writes their drafts;
// `options.readOtherModeArtifact` returns the other mode's artifact slot that
// the root keeps (E1), `options.beforeExport` learns whether Save projects
// the draft request, and `options.onSaved` runs once the download is handed off.
/**
 * What the composition root passes to `exportSession`.
 * @typedef {object} SessionSaveOptions
 * @property {(operation: 'save' | 'load') => SessionBusy | null} [availability]
 * @property {DownloadChoicePort} [askChoice]
 * @property {{ value?: Parameters<typeof bindRecordDisplayDrafts>[1] }} [recordDisplayRows]
 * @property {() => Readonly<ArtifactSlot> | null} [readOtherModeArtifact]
 * @property {() => Promise<string | null>} [resolveTitle] Null when the person canceled the title dialog.
 * @property {(request: { draftRequest: boolean }) => Promise<Pick<ExportSessionDocumentOptions,
 *   'linearRecordCatalog' | 'omittedAnnotationRecordKeys'>>} [beforeExport]
 * @property {(result: Extract<SessionSaveResult, { status: 'saved' }>) => void} [onSaved]
 * @property {(error: UserFacingError | null) => void} [onError]
 */

/**
 * @param {string | null} [titleOverride]
 * @param {SessionSaveOptions} [options]
 */
export const exportSession = (titleOverride = null, options = {}) => {
  if (sessionSaveInFlight) {
    recordSessionLifecycleEvent('session-save-joined');
    return sessionSaveInFlight.promise;
  }
  const { availability = sessionOperationAvailability } = options;
  const busy = availability('save');
  if (busy) return Promise.resolve(busy);
  /** @type {{ canceled: boolean, promise: Promise<any> | null }} */
  const operation = { canceled: false, promise: null };
  const previousAlert = state.errorLog.value;
  const isCurrent = () => sessionSaveInFlight === operation && !operation.canceled;
  operation.promise = Promise.resolve().then(async () => {
    const busy = availability('save');
    if (busy) return busy;
    if (!isCurrent()) return { status: 'canceled' };
    const title = options.resolveTitle ? await options.resolveTitle() : titleOverride;
    if (title === null && options.resolveTitle) return;
    if (!isCurrent()) return { status: 'canceled' };
    state.sessionSavePending.value = true;
    recordSessionLifecycleEvent('session-save-pending-published');
    // Only mutable draft configuration/navigation is copied. Adopted biological
    // payloads, resources, catalogs, caches and Results retain their existing owner.
    // Session 46 keeps each drawing as its mode's slice (PD-OI-086); the
    // Features-list record is the shown Result's. `drawing` is the shown mode's.
    // The cast: outside a Load the shown mode is a drawing's (see `captureSessionImportSnapshot`).
    const drawing = state.drawings[/** @type {DiagramMode} */ (state.mode.value)];
    // Each drawing's drafts are written on the rows of their records (OV-278,
    // OV-325): a Session without a Web draft loads a draft selected by record ID on '#1'.
    const recordDisplayRows = /** @type {{ recordDisplayRows?: { value?: Parameters<typeof bindRecordDisplayDrafts>[1] } }} */ (options).recordDisplayRows?.value || [];
    /** @param {'circular' | 'linear'} mode @returns {DrawingState} */
    const savedDrawing = (mode) => {
      const live = state.drawings[mode];
      const drafts = bindRecordDisplayDrafts(live.recordDisplayDrafts, recordDisplayRows);
      return drafts === live.recordDisplayDrafts ? live : { ...live, recordDisplayDrafts: drafts };
    };
    // Each drawing is checked as it is, before the JSON copy, so a value JSON
    // cannot write (Infinity, NaN) fails Save instead of being written as null.
    SLICE_MODES.forEach((mode) => validateCurrentWriterActiveConfig({
      mode, storedConfig: buildConfigData(savedDrawing(mode))
    }));
    /** @param {'circular' | 'linear'} mode */
    const sliceOf = (mode) => buildModeSliceData(savedDrawing(mode), mode, {
      selectedFeatureRecordIdx: mode === state.mode.value ? Number(state.selectedFeatureRecordIdx.value) || 0 : 0
    });
    const modes = { circular: sliceOf('circular'), linear: sliceOf('linear') };
    SLICE_MODES.forEach((mode) => validateModeSliceFields(modes[mode]));
    const readArtifacts = () => chooseSessionArtifacts(
      /** @type {{ readOtherModeArtifact?: () => Readonly<ArtifactSlot> | null }} */ (options).readOtherModeArtifact?.() ?? null
    );
    let { artifact, otherArtifact } = readArtifacts();
    if (!artifact.committedCanonicalSession
      && hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs })) {
      assertActiveModeInputs();
    }
    const artifactUi = artifactUiState(artifact);
    const savedUi = {
      mode: state.mode.value,
      ...previewViewportToPersist(),
      selectedResultIndex: artifactUi.selectedResultIndex,
      generatedLegendPosition: artifactUi.generatedLegendPosition,
      generatedMultiRecordCanvas: artifactUi.generatedMultiRecordCanvas,
      generatedCircularPlotTitlePosition: artifactUi.generatedCircularPlotTitlePosition,
      featurePanelTab: state.featurePanelTab.value,
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      downloadDpi: state.downloadDpi.value,
      autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
      paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
      appliedPaletteName: artifactUi.appliedPaletteName,
      appliedPaletteColors: artifactUi.appliedPaletteColors,
      // App-level settings shared by both drawings.
      losatExecution: buildLosatExecutionData(),
      richFeaturePopup: Boolean(state.richFeaturePopup.value)
    };
    const prepared = await options.beforeExport?.({ draftRequest: !artifact.committedCanonicalSession });
    if (!isCurrent()) return { status: 'canceled' };
    // `beforeExport` may show edits on the displayed Result, which commits
    // its content: the Results are read as they are now.
    ({ artifact, otherArtifact } = readArtifacts());
    const result = await exportSessionDocument(title, {
      ...options, ...prepared, drawing, modes, savedUi, isCurrent, artifact, otherArtifact
    });
    if (result?.status === 'saved') options.onSaved?.(result);
    return result;
  }).catch((error) => {
    recordSessionLifecycleEvent('session-save-error');
    if (!isCurrent() || (state.errorLog.value !== previousAlert && state.errorLog.value !== null)) {
      return { status: 'stale' };
    }
    if (typeof options.onError !== 'function') throw error;
    const model = normalizeUserFacingError(error, { operation: 'session-save' });
    options.onError(model);
    return { status: 'error', error: model };
  }).finally(() => {
    if (sessionSaveInFlight === operation) {
      state.sessionSavePending.value = false;
      sessionSaveInFlight = null;
    }
    recordSessionLifecycleEvent('session-save-pending-cleared');
  });
  sessionSaveInFlight = operation;
  return operation.promise;
};

/**
 * The Session file input's change event, or an object with its `target`.
 * @typedef {{ target?: { files?: ArrayLike<File> | null, value: string } | null } | null | undefined} SessionFileInput
 */
/**
 * What the composition root passes to `importSession`: the Save and Load availability (as `SessionSaveOptions`
 * has it) and a step to run before the Load starts, both read here; the required Load port
 * `confirmLoadRepairs`, which shows the list of what the file cannot fully restore and answers [Load] or
 * [Cancel] (every caller that loads decides; there is no default); every other option is passed on to
 * `importSessionDocument` unread.
 * @typedef {Pick<SessionSaveOptions, 'availability'> & { beforeImport?: () => void | Promise<void>,
 *   confirmLoadRepairs: ConfirmLoadRepairs } & Record<string, unknown>} SessionLoadOptions
 */
/**
 * @param {SessionFileInput} event
 * @param {SessionLoadOptions} options
 */
export const importSession = async (event, options) => {
  const input = event?.target;
  const file = input?.files?.[0];
  if (!file) return { status: 'skipped' };
  // A caller without the port loads a file Load restores fully; a file with
  // repairs fails its Load (the port is required).
  const { availability = sessionOperationAvailability } = options ?? {};
  const busy = availability('load');
  if (busy) {
    input.value = '';
    return busy;
  }
  const operation = { canceled: false, controller: new AbortController() };
  activeSessionImport = operation;
  const isCurrent = () => activeSessionImport === operation && !operation.canceled;
  state.sessionImportPending.value = true;
  recordSessionLifecycleEvent('session-import-pending-published');
  try {
    await options?.beforeImport?.();
    if (!isCurrent()) return { status: 'canceled' };
    const result = await importSessionDocument({ target: { files: [file], value: input.value } }, {
      ...options, isCurrent, signal: operation.controller.signal
    });
    return result;
  } finally {
    if (activeSessionImport === operation) {
      activeSessionImport = null;
      state.sessionImportPending.value = false;
    }
    input.value = '';
    recordSessionLifecycleEvent('session-import-pending-cleared');
    // Discovery watchers see the cleared pending flag while an adopted
    // current-schema Session still defers resource reads; clear it afterwards.
    await nextTick();
    if (!activeSessionImport) state.sessionResourceDiscoveryDeferred.value = false;
  }
};
