// @ts-check
import { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat, validateCurrentWriterActiveConfig } from './session-active-config-contract.js';
import { resolveActiveLayoutPreference } from '../app/layout-preferences.js';
import { migrateLegacyLinearLabelVisibility } from '../app/linear-label-visibility.js';
import { migrateLegacyRecordDisplayDrafts } from '../app/record-display-options.js';
import { FEATURE_CATALOG_SCHEMA, migrateLegacyFeatureCatalog } from './feature-catalog.js';
import { migrateSessionFeatureEdits, migrateSessionFeaturePlacements } from './feature-edit-migration.js';
import { adoptCurrentSessionResources } from './session-resource-backing.js';
import { defaultFeatureRendering } from '../utils/feature-rendering.js';
const CURRENT_VERSION = 45, CURRENT_REQUEST_SCHEMA = 9, ACCEPTED_REQUEST_SCHEMAS = new Set([CURRENT_REQUEST_SCHEMA]), HISTORICAL_VERSIONS = new Set([31, 32, 33, 39]), CACHE_LIMIT_BYTES = 64 * 1024 * 1024;
const ARTIFACT_FIELDS = ['results', 'features', 'editorState', 'orthogroupState', 'runMetadata', 'losatCache', 'losatDerivedCache', 'proteinIdentityManifest'];
// A published Gallery file carries no draft intent for its unused mode (GUI
// remediation S00 decision 1). These fields, read only by the other mode, are
// written with fresh defaults; mode profiles own the per-mode title and fonts.
const UNUSED_MODE_FRESH_FIELDS = Object.freeze({
  circular: Object.freeze({ adv: Object.freeze(['linear_show_replicon', 'linear_accession_visibility', 'linear_length_visibility']) }),
  linear: Object.freeze({ form: Object.freeze(['multi_record_canvas']) })
});
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
const isCliWritten = (session) => !has(session, 'config') && session?.cliInvocation?.generatedBy === 'gbdraw';
const validateEnvelope = (session) => {
  if (!isObject(session) || session.format !== 'gbdraw-session') throw new Error('Gallery publication requires a gbdraw-session document.'); if (Number(session.version) !== CURRENT_VERSION) throw new Error(`Gallery publication requires session version ${CURRENT_VERSION}.`);
  if (!isObject(session.renderRequest) || !ACCEPTED_REQUEST_SCHEMAS.has(Number(session.renderRequest.schema))) throw new Error(`Gallery publication requires canonical renderRequest schema ${CURRENT_REQUEST_SCHEMA}.`);
  return session;
};
const validateCurrent = (session) => {
  validateEnvelope(session); validateCurrentWriterActiveConfig({ mode: session.renderRequest.mode, storedConfig: session.config });
  return session;
};
const publicationConfig = (session, projection) => {
  const projected = clone(projection.config || {}), stored = session.config || {}, losat = createDefaultLosat();
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
  for (const key of ['palette', 'annotationSets', 'recordDisplayDrafts', 'featurePlacementOverrides', 'modeProfiles', 'linearRecordLayout', 'losatProgram'])
    if (has(stored, key)) config[key] = clone(stored[key]);
  if (stored.colorsAreOverrides === true && Object.keys(stored.colors || {}).length) Object.assign(config, { colors: clone(stored.colors), colorsAreOverrides: true });
  const colors = session.renderRequest.diagramOptions?.colors;
  if (!colors?.defaultColors && !colors?.defaultColorsFile) Object.assign(config, { colors: {}, colorsAreOverrides: false });
  for (const key of ['webEdits', 'paletteInstantPreviewEnabled']) if (has(stored, key)) config[key] = clone(stored[key]);
  delete config.blastSource; delete config.adv.losatProgram;
  if (isCliWritten(session)) {
    // Legend and plot-title positions are layout preferences of the request.
    const committedLayout = resolveActiveLayoutPreference(projection.layoutPreferences, projection.mode,
      Boolean(config.form.multi_record_canvas));
    Object.assign(config.form, { legend: committedLayout.legend });
    Object.assign(config.adv, { plot_title_position: committedLayout.plotTitlePosition });
  }
  const fresh = { form: createDefaultForm(), adv: createDefaultAdv(projection.mode) };
  for (const [domain, fields] of Object.entries(UNUSED_MODE_FRESH_FIELDS[projection.mode]))
    for (const field of fields) config[domain][field] = fresh[domain][field];
  return config;
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
    resources: session.resources, webFiles: session.webFiles || {}, legacyFiles: session.files, storedConfig: session.config,
    initializeCliInputs: isCliWritten(session),
    fileBindings: session.cliInvocation?.fileBindings, sessionResourceTable: adoptCurrentSessionResources(session.resources),
    deferResourceContent: false, adoptCanonicalPayloads: true });
  const config = publicationConfig(session, projection); validateCurrentWriterActiveConfig({ mode: projection.mode, storedConfig: config });
  const filesData = projection.files;
  if (projection.mode === 'linear') filesData.linearSeqs.forEach((sequence, index) => {
    sequence.cardinality = renderRequest.records[index]?.cardinality;
  });
  const state = owners.buildRequestState({ session, projection, config, filesData });
  const plan = projection.mode === 'linear' ? owners.resolveComparisonPlan({ plan: state.linearComparisonPlan, sequences: filesData.linearSeqs,
    layout: state.linearRecordLayoutEnabled.value ? state.linearRecordRows : [],
    losatProgram: state.losatProgram.value, blastpMode: state.losat?.blastp?.mode }) : null;
  const rebuilt = owners.buildRequest({ state, filesData, comparisonPlanSnapshot: plan });
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
  return { config, rebuilt, equivalence: await owners.assertRequestsEquivalent({ expectedRequest: comparable(renderRequest),
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
 * @property {(input: Record<string, any>) => Record<string, any>} buildRequestState
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
    // Session 45 keys per-feature edits by source identity (design Q4 4.3).
    features: migrateSessionFeatureEdits({
      features: session?.features,
      mode: session?.renderRequest?.mode,
      catalog: session?.editorState?.featureCatalog || null
    }).features
  });
  const admit = (session) => {
    const version = Number(session?.version);
    if (version === CURRENT_VERSION) {
      if ([7, 8].includes(Number(session.renderRequest?.schema))) {
        return validateCurrent({ ...session, renderRequest: publicationCanonicalRequest(
          session.renderRequest,
          owners.promoteRequest,
          {
            featureCatalog: session.editorState?.featureCatalog || null,
            legacyOrthogroupState: session.orthogroupState || null
          }
        ) });
      }
      return isCliWritten(session) ? validateEnvelope(session) : validateCurrent(session);
    }
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
    const admitted = admit(source), { config, rebuilt, equivalence } = await rebuild(admitted);
    const session = {
      ...admitted,
      config,
      renderRequest: rebuilt.renderRequest,
      resources: rebuilt.resources,
      webFiles: rebuilt.webFiles,
      orthogroupState: currentOrthogroupState(admitted.orthogroupState)
    };
    return { session: validateCurrent(session), equivalence };
  };
  const finalize = async ({ prepared, replayed }) => {
    validateCurrent(prepared); validateCurrent(replayed);
    const resources = mergeReplayResources(prepared, replayed);
    await owners.assertRequestsEquivalent({ expectedRequest: prepared.renderRequest, expectedResources: prepared.resources, actualRequest: replayed.renderRequest,
      actualResources: replayed.resources, normalizeReplayGeneratedResources: true });
    let session = { ...prepared, resources };
    for (const field of ARTIFACT_FIELDS) if (has(replayed, field)) session[field] = replayed[field];
    session = applyDerivedCachePublicationPolicy(session);
    return { session, equivalence: (await rebuild(session)).equivalence };
  };
  const validate = (session) => (validateCurrent(session), rebuild(session));
  return Object.freeze({ admitGallerySession: admit, finalizeGallerySessionPublication: finalize, prepareGallerySessionForPublication: prepare, validateGalleryPublicationReadiness: validate });
};
