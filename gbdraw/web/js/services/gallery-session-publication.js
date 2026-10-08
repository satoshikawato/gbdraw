// @ts-check
import { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat, validateCurrentWriterActiveConfig } from './session-active-config-contract.js';
import {
  LEGACY_LAYOUT_PREFERENCE_FIELDS, normalizeLayoutPreferences, resolveActiveLayoutPreference, restoredLayoutPreferences
} from './layout-preferences.js';
import { splitDraftIntoModes, validateModeSliceFields } from './mode-scoped-migration.js';
import { migrateLegacyLinearLabelVisibility } from './linear-label-visibility.js';
import { migrateLegacyRecordDisplayDrafts } from './record-display-model.js';
import { FEATURE_CATALOG_SCHEMA, migrateLegacyFeatureCatalog } from './feature-catalog.js';
import { migrateSessionFeatureEdits, migrateSessionFeaturePlacements } from './feature-edit-migration.js';
import { adoptCurrentSessionResources } from './session-resource-backing.js';
import { defaultFeatureRendering } from '../utils/feature-rendering.js';
const CURRENT_VERSION = 46, CURRENT_REQUEST_SCHEMA = 9, ACCEPTED_REQUEST_SCHEMAS = new Set([CURRENT_REQUEST_SCHEMA]), HISTORICAL_VERSIONS = new Set([31, 32, 33, 39]), CACHE_LIMIT_BYTES = 64 * 1024 * 1024;
const ARTIFACT_FIELDS = ['results', 'editorState', 'orthogroupState', 'runMetadata', 'losatCache', 'losatDerivedCache', 'proteinIdentityManifest'];
// The top-level homes of the flat draft that Session 46 moved into `modes`.
const RETIRED_DRAFT_FIELDS = ['config', 'features'];
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
const clone = (value) => value === undefined ? undefined : JSON.parse(JSON.stringify(value));
const has = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
const currentOrthogroupState = (value) => {
  if (!isObject(value)) return value;
  const current = clone(value);
  delete current.selectedOrthogroupAlignmentFeature;
  return current;
};
const regenerableProteinCache = (session) => session.renderRequest?.comparisons?.some(
  ({ kind }) => kind === 'generatedProteinComparison') && session.proteinIdentityManifest?.schema === 2
  && session.losatCache?.entries?.some((entry) => entry?.schema === 4 && entry?.kind === 'raw-losat' && entry?.program === 'blastp'
    && entry?.idEncoding === 'runtime-handle-v1'
    && typeof entry?.queryProteinSetHash === 'string' && typeof entry?.subjectProteinSetHash === 'string');
/**
 * @param {Record<string, any>} session
 * @param {{ limitBytes?: number }} [options]
 * @returns {Record<string, any>}
 */
export const applyDerivedCachePublicationPolicy = (session, { limitBytes = CACHE_LIMIT_BYTES } = {}) => {
  const entries = session.losatDerivedCache?.entries;
  if (!entries?.length || !regenerableProteinCache(session) || new TextEncoder().encode(JSON.stringify(entries)).byteLength <= limitBytes) return session;
  return { ...session, losatDerivedCache: { ...session.losatDerivedCache, entries: [] } };
};
// A CLI-written Session has no Web configuration; publication derives it from
// the request, as Session Load does, so a Gallery Session can be built from its
// declared command.
const isCliWritten = (session) => !has(session, 'config') && !has(session, 'modes')
  && session?.cliInvocation?.generatedBy === 'gbdraw';
// The draft of the Session's mode: its slice (Session 46) or the flat draft of
// an older Session that admission promoted. Publication writes the slice of the
// Session's mode only; the other mode takes its defaults (plan 4.3).
const modeSlice = (session) => {
  const slice = isObject(session?.modes) ? session.modes[session.renderRequest?.mode] : null;
  return isObject(slice) ? slice : null;
};
const draftConfig = (session) => (has(session, 'modes') ? modeSlice(session)?.config || {} : session?.config);
const draftFeatures = (session) => (has(session, 'modes') ? modeSlice(session)?.features : session?.features);
const validateEnvelope = (session) => {
  if (!isObject(session) || session.format !== 'gbdraw-session') throw new Error('Gallery publication requires a gbdraw-session document.'); if (Number(session.version) !== CURRENT_VERSION) throw new Error(`Gallery publication requires session version ${CURRENT_VERSION}.`);
  if (!isObject(session.renderRequest) || !ACCEPTED_REQUEST_SCHEMAS.has(Number(session.renderRequest.schema))) throw new Error(`Gallery publication requires canonical renderRequest schema ${CURRENT_REQUEST_SCHEMA}.`);
  // A Gallery Session shows one diagram mode (E1: `otherModeResult` holds a second).
  if (has(session, 'otherModeResult')) throw new Error('Gallery publication requires a Session with one diagram mode.');
  return session;
};
// A Session 46 keeps the Web draft in `modes` only (plan 4.1).
const validateCurrentFields = (session) => {
  if (RETIRED_DRAFT_FIELDS.some((field) => has(session, field))) throw new Error('Gallery publication requires the Web draft in modes (Session 46).');
  if (has(session, 'modes')) {
    const modes = session.modes;
    if (!isObject(modes) || Object.keys(modes).some((mode) => mode !== 'circular' && mode !== 'linear'))
      throw new Error('Gallery publication requires modes to hold a circular and a linear slice only.');
    Object.values(modes).forEach((slice) => validateModeSliceFields(slice));
  }
  return session;
};
const validateCurrent = (session) => {
  validateEnvelope(session); validateCurrentWriterActiveConfig({ mode: session.renderRequest.mode, storedConfig: draftConfig(session) });
  return session;
};
const publicationConfig = (session, projection) => {
  const projected = clone(projection.config || {}), stored = draftConfig(session) || {}, losat = createDefaultLosat();
  const config = { ...projected,
    form: { ...createDefaultForm(), ...projected.form, ...stored.form },
    adv: { ...createDefaultAdv(projection.mode), ...projected.adv, ...stored.adv },
    losat: { ...losat, ...projected.losat, ...stored.losat,
      blastn: { ...losat.blastn, ...projected.losat?.blastn, ...stored.losat?.blastn },
      blastp: { ...losat.blastp, ...projected.losat?.blastp, ...stored.losat?.blastp } },
    circularConservation: { ...createDefaultCircularConservation(), ...projected.circularConservation, ...clone(stored.circularConservation) },
    linearComparisonPlan: clone(stored.linearComparisonPlan || projected.linearComparisonPlan
      || { mode: 'none', defaultSource: 'losat', edges: [] })
  };
  for (const key of ['palette', 'annotationSets', 'recordDisplayDrafts', 'featurePlacementOverrides', 'linearRecordLayout', 'losatProgram'])
    if (has(stored, key)) config[key] = clone(stored[key]);
  if (stored.colorsAreOverrides === true && Object.keys(stored.colors || {}).length) Object.assign(config, { colors: clone(stored.colors), colorsAreOverrides: true });
  const colors = session.renderRequest.diagramOptions?.colors;
  if (!colors?.defaultColors && !colors?.defaultColorsFile) Object.assign(config, { colors: {}, colorsAreOverrides: false });
  for (const key of ['webEdits', 'paletteInstantPreviewEnabled']) if (has(stored, key)) config[key] = clone(stored[key]);
  delete config.blastSource; delete config.adv.losatProgram;
  return config;
};
// The layout preferences of the Session's drawing: a Session 46 slice's slot
// over its request's layout, or an older Session's layout as Session Load
// restores it. The request reads the shown slot's Legend and plot-title
// positions from the drawing's form and advanced settings.
const publicationLayout = (session, config, projection) => {
  const mode = projection.mode;
  const slot = modeSlice(session)?.ui?.layoutPreferences;
  const ui = isObject(session.ui) ? session.ui : {};
  const layout = has(session, 'modes')
    ? normalizeLayoutPreferences({ ...normalizeLayoutPreferences(projection.layoutPreferences),
      ...(isObject(slot) ? { [mode]: slot } : {}) })
    : restoredLayoutPreferences(ui, { mode, multiRecord: Boolean(config.form?.multi_record_canvas),
      projected: projection.layoutPreferences,
      active: { legend: config.form?.legend, plotTitlePosition: config.adv?.plot_title_position } });
  const active = resolveActiveLayoutPreference(layout, mode, Boolean(config.form?.multi_record_canvas));
  Object.assign(config.form, { legend: active.legend });
  Object.assign(config.adv, { plot_title_position: active.plotTitlePosition });
  return layout;
};
// The published Session 46: the Session's mode's slice holds the rebuilt
// configuration and its layout slot; the other slice is absent, so it takes
// that mode's defaults. An older Session's flat draft is split by the registry
// (plan 4.2); a Session 46 keeps the rest of its slice. The layout has one
// home, so the older `ui` layout fields go.
const publishedSession = (session, config, layout, projection) => {
  const mode = projection.mode;
  const flat = !has(session, 'modes');
  const draft = flat ? { ...session, config, ui: { ...(isObject(session.ui) ? session.ui : {}), mode } } : { config, ui: { mode } };
  // Override colors merge into the committed request's colors.
  const split = splitDraftIntoModes(draft, { committedMode: mode, modeProfiles: null, paletteColors: projection.config?.colors || null });
  const slice = flat ? split.modes[mode] : { ...modeSlice(session), config: split.modes[mode].config };
  const published = flat ? split : { ...session };
  const ui = Object.fromEntries(Object.entries(isObject(published.ui) ? published.ui : {})
    .filter(([field]) => !LEGACY_LAYOUT_PREFERENCE_FIELDS.includes(field)));
  return { ...published, ui, modes: { [mode]: { ...slice, ui: { ...(slice.ui || {}), layoutPreferences: clone(layout[mode]) } } } };
};
const publicationCanonicalRequest = (
  request,
  promoteRequest,
  { featureCatalog = null, legacyOrthogroupState = null } = {}
) => {
  const current = promoteRequest(request, { featureCatalog, legacyOrthogroupState });
  if (Number(request?.schema) === 5) {
    current.records.forEach((record) => { record.cardinality = 'exactly_one'; });
  }
  return current;
};
const rebuildIntent = async (session, owners) => {
  const renderRequest = publicationCanonicalRequest(
    session.renderRequest,
    owners.promoteRequest,
    {
      featureCatalog: session.editorState?.featureCatalog || null,
      legacyOrthogroupState: session.orthogroupState || null
    }
  );
  const projection = owners.projectRequest({ renderRequest,
    resources: session.resources, webFiles: session.webFiles || {}, legacyFiles: session.files, storedConfig: draftConfig(session),
    initializeCliInputs: isCliWritten(session),
    fileBindings: session.cliInvocation?.fileBindings, sessionResourceTable: adoptCurrentSessionResources(session.resources),
    deferResourceContent: false, adoptCanonicalPayloads: true });
  const config = publicationConfig(session, projection); validateCurrentWriterActiveConfig({ mode: projection.mode, storedConfig: config });
  const layout = publicationLayout(session, config, projection);
  const filesData = projection.files;
  if (projection.mode === 'linear') filesData.linearSeqs.forEach((sequence, index) => {
    sequence.cardinality = renderRequest.records[index]?.cardinality;
  });
  const { state, drawing } = owners.buildRequestState({ session: { ...session, features: draftFeatures(session) }, projection, config, filesData });
  const plan = projection.mode === 'linear' ? owners.resolveComparisonPlan({ plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs,
    layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
    losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat?.blastp?.mode }) : null;
  const rebuilt = owners.buildRequest({ state, drawing, filesData, comparisonPlanSnapshot: plan });
  if (!isObject(rebuilt.renderRequest.output) || !isObject(session.renderRequest.output)) throw new Error('Gallery publication cannot preserve committed output metadata policy.');
  rebuilt.renderRequest.output.interactiveMetadataPolicy = session.renderRequest.output.interactiveMetadataPolicy;
  // A CLI-written request carries the resolved configuration; publication
  // writes the Web's own configOverrides, so Session Load needs no Worker. The
  // refresh tool checks that the replayed figure equals the declared figure.
  const cliConfig = isCliWritten(session) && isObject(session.renderRequest.diagramOptions?.config);
  if (isObject(session.renderRequest.diagramOptions?.config) && !cliConfig) {
    rebuilt.renderRequest.diagramOptions.config = clone(session.renderRequest.diagramOptions.config); delete rebuilt.renderRequest.diagramOptions.configOverrides;
  }
  // The CLI omits feature renderings it leaves at their defaults; the Web
  // request states them (repeat_region underlay), so defaults compare as absent.
  const comparable = (request) => {
    if (!cliConfig) return request;
    const diagramOptions = Object.fromEntries(Object.entries(request.diagramOptions)
      .filter(([key]) => key !== 'config' && key !== 'configOverrides'));
    const featureShapes = Object.fromEntries(Object.entries(diagramOptions.featureShapes || {})
      .filter(([type, rendering]) => rendering !== defaultFeatureRendering(type)));
    if (Object.keys(featureShapes).length) diagramOptions.featureShapes = featureShapes;
    else delete diagramOptions.featureShapes;
    return { ...request, diagramOptions };
  };
  return { config, layout, projection, rebuilt, equivalence: await owners.assertRequestsEquivalent({ expectedRequest: comparable(renderRequest),
    expectedResources: session.resources, actualRequest: comparable(rebuilt.renderRequest), actualResources: rebuilt.resources }) };
};
const mergeReplayResources = (prepared, replayed) => {
  for (const [id, expected] of Object.entries(prepared.resources || {})) {
    const actual = replayed.resources?.[id], fields = actual && ['kind', 'encoding', 'data'].filter((field) => expected[field] !== actual[field]);
    if (fields?.length) throw new Error(`Gallery replay resource collision at resources.${id}: ${fields.join(', ')}.`);
  }
  const resources = { ...(prepared.resources || {}) }, referenced = new Set();
  const collect = (value) => {
    if (Array.isArray(value)) return value.forEach(collect);
    if (!isObject(value)) return;
    if (typeof value.resourceId === 'string' && value.resourceId) referenced.add(value.resourceId); Object.values(value).forEach(collect);
  };
  ARTIFACT_FIELDS.forEach((field) => collect(replayed[field]));
  for (const id of referenced) {
    if (resources[id]) continue;
    if (!replayed.resources?.[id]) throw new Error(`Gallery replay artifact references missing resource '${id}'.`);
    resources[id] = replayed.resources[id];
  }
  return resources;
};
/**
 * The owners a publication calls: the Session and request owners
 * (`services/session-request.js`, `services/gallery-session-migration.js`) and
 * the Linear comparison planner. Each takes and returns unvalidated documents.
 * @typedef {object} GallerySessionPublicationOwners
 * @property {(session: Record<string, any>) => Record<string, any>} promoteSession
 *   Promotes a historical Session to the current version.
 * @property {(input: Record<string, any>) => Promise<Record<string, any>>} assertRequestsEquivalent
 * @property {(input: Record<string, any>) => Record<string, any>} buildRequest
 * @property {(input: Record<string, any>) => { state: Record<string, any>, drawing: Record<string, any> }} buildRequestState
 *   The request inputs of a Session: project inputs and artifacts (`state`) and settings and edits (`drawing`).
 * @property {(request: Record<string, any>, promotion?: Record<string, any>) => Record<string, any>} promoteRequest
 * @property {(input: Record<string, any>) => Record<string, any>} projectRequest
 * @property {(input: Record<string, any>) => any} resolveComparisonPlan
 */

/**
 * @param {GallerySessionPublicationOwners} owners
 */
export const createGallerySessionPublication = (owners) => {
  const promoteVisibilityState = (session) => ({
    ...session,
    config: isObject(session?.config) && isObject(session.config.adv)
      ? {
          ...session.config,
          adv: migrateLegacyLinearLabelVisibility(session.config.adv),
          ...(Array.isArray(session.config.recordDisplayDrafts) ? {
            recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(
              session.config.recordDisplayDrafts
            )
          } : {}),
          ...(has(session.config, 'featurePlacementOverrides') ? {
            featurePlacementOverrides: migrateSessionFeaturePlacements(session.config.featurePlacementOverrides)
          } : {})
        }
      : session?.config,
    editorState: [3, 4].includes(session?.editorState?.featureCatalog?.schema)
      && session.editorState.featureCatalog.schema !== FEATURE_CATALOG_SCHEMA
      ? {
          ...session.editorState,
          featureCatalog: migrateLegacyFeatureCatalog(
            session.editorState.featureCatalog
          )
        }
      : session?.editorState,
    // The current writer keys per-feature edits by source identity (design Q4 4.3).
    features: migrateSessionFeatureEdits({
      features: session?.features,
      mode: session?.renderRequest?.mode,
      catalog: session?.editorState?.featureCatalog || null
    }).features
  });
  const admit = (session) => {
    const version = Number(session?.version);
    if (version === CURRENT_VERSION) {
      validateCurrentFields(session);
      return isCliWritten(session) ? validateEnvelope(session) : validateCurrent(session);
    }
    // An older Session is promoted to the current request and kept as a flat
    // draft until `prepare` splits it into the Session's mode slice.
    if ([40, 41, 42, 44].includes(version)) return validateCurrent(promoteVisibilityState({ ...session, version: CURRENT_VERSION,
      renderRequest: publicationCanonicalRequest(
        session.renderRequest,
        owners.promoteRequest,
        {
          featureCatalog: session.editorState?.featureCatalog || null,
          legacyOrthogroupState: session.orthogroupState || null
        }
      ) }));
    if (!HISTORICAL_VERSIONS.has(version)) throw new Error(`Gallery publication supports current version ${CURRENT_VERSION} or historical versions 31-33/39-44; received ${String(session?.version)}.`);
    return validateCurrent(promoteVisibilityState(owners.promoteSession(session)));
  };
  const rebuild = (session) => rebuildIntent(session, owners);
  const prepare = async (source) => {
    const admitted = admit(source), { config, layout, projection, rebuilt, equivalence } = await rebuild(admitted);
    const session = {
      ...publishedSession(admitted, config, layout, projection),
      renderRequest: rebuilt.renderRequest,
      resources: rebuilt.resources,
      webFiles: rebuilt.webFiles,
      orthogroupState: currentOrthogroupState(admitted.orthogroupState)
    };
    return { session: validateCurrent(validateCurrentFields(session)), equivalence };
  };
  const finalize = async ({ prepared, replayed }) => {
    validateCurrent(validateCurrentFields(prepared)); validateCurrent(validateCurrentFields(replayed));
    const resources = mergeReplayResources(prepared, replayed);
    await owners.assertRequestsEquivalent({ expectedRequest: prepared.renderRequest, expectedResources: prepared.resources, actualRequest: replayed.renderRequest,
      actualResources: replayed.resources, normalizeReplayGeneratedResources: true });
    let session = { ...prepared, resources };
    for (const field of ARTIFACT_FIELDS) if (has(replayed, field)) session[field] = replayed[field];
    session = applyDerivedCachePublicationPolicy(session);
    return { session, equivalence: (await rebuild(session)).equivalence };
  };
  const validate = (session) => (validateCurrent(validateCurrentFields(session)), rebuild(session));
  return Object.freeze({ admitGallerySession: admit, finalizeGallerySessionPublication: finalize, prepareGallerySessionForPublication: prepare, validateGalleryPublicationReadiness: validate });
};
