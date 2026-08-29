import { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat, validateCurrentWriterActiveConfig } from './session-active-config-contract.js';
import { adoptCurrentSessionResources } from './session-resource-backing.js';
const CURRENT_VERSION = 40, CURRENT_REQUEST_SCHEMA = 6, ACCEPTED_REQUEST_SCHEMAS = new Set([5, CURRENT_REQUEST_SCHEMA]), HISTORICAL_VERSIONS = new Set([31, 32, 33, 39]), CACHE_LIMIT_BYTES = 64 * 1024 * 1024;
const ARTIFACT_FIELDS = ['results', 'features', 'editorState', 'orthogroupState', 'runMetadata', 'losatCache', 'losatDerivedCache', 'proteinIdentityManifest'];
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
const clone = (value) => value === undefined ? undefined : JSON.parse(JSON.stringify(value));
const has = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
const regenerableProteinCache = (session) => session.renderRequest?.comparisons?.some(
  ({ kind }) => kind === 'generatedProteinComparison') && session.proteinIdentityManifest?.schema === 2
  && session.losatCache?.entries?.some((entry) => entry?.schema === 4 && entry?.kind === 'raw-losat' && entry?.program === 'blastp'
    && entry?.idEncoding === 'runtime-handle-v1'
    && typeof entry?.queryProteinSetHash === 'string' && typeof entry?.subjectProteinSetHash === 'string');
export const applyDerivedCachePublicationPolicy = (session, { limitBytes = CACHE_LIMIT_BYTES } = {}) => {
  const entries = session.losatDerivedCache?.entries;
  if (!entries?.length || !regenerableProteinCache(session) || new TextEncoder().encode(JSON.stringify(entries)).byteLength <= limitBytes) return session;
  return { ...session, losatDerivedCache: { ...session.losatDerivedCache, entries: [] } };
};
const validateCurrent = (session) => {
  if (!isObject(session) || session.format !== 'gbdraw-session') throw new Error('Gallery publication requires a gbdraw-session document.'); if (Number(session.version) !== CURRENT_VERSION) throw new Error(`Gallery publication requires session version ${CURRENT_VERSION}.`);
  if (!isObject(session.renderRequest) || !ACCEPTED_REQUEST_SCHEMAS.has(Number(session.renderRequest.schema))) throw new Error('Gallery publication requires canonical renderRequest schema 5 or 6.'); validateCurrentWriterActiveConfig({ mode: session.renderRequest.mode, storedConfig: session.config });
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
  for (const key of ['palette', 'annotationSets', 'modeProfiles', 'linearRecordLayout', 'losatProgram'])
    if (has(stored, key)) config[key] = clone(stored[key]);
  if (stored.colorsAreOverrides === true && Object.keys(stored.colors || {}).length) Object.assign(config, { colors: clone(stored.colors), colorsAreOverrides: true });
  const colors = session.renderRequest.diagramOptions?.colors;
  if (!colors?.defaultColors && !colors?.defaultColorsFile) Object.assign(config, { colors: {}, colorsAreOverrides: false });
  for (const key of ['webEdits', 'paletteInstantPreviewEnabled']) if (has(stored, key)) config[key] = clone(stored[key]);
  delete config.blastSource; delete config.adv.losatProgram;
  return config;
};
const publicationCanonicalRequest = (request, promoteRequest) => {
  const current = promoteRequest(request);
  if (Number(request?.schema) === 5) {
    current.records.forEach((record) => { record.cardinality = 'exactly_one'; });
  }
  return current;
};
const rebuildIntent = async (session, owners) => {
  const renderRequest = publicationCanonicalRequest(session.renderRequest, owners.promoteRequest);
  const projection = owners.projectRequest({ renderRequest,
    resources: session.resources, webFiles: session.webFiles || {}, legacyFiles: session.files, storedConfig: session.config,
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
  if (isObject(session.renderRequest.diagramOptions?.config)) {
    rebuilt.renderRequest.diagramOptions.config = clone(session.renderRequest.diagramOptions.config); delete rebuilt.renderRequest.diagramOptions.configOverrides;
  }
  return { config, rebuilt, equivalence: await owners.assertRequestsEquivalent({ expectedRequest: renderRequest,
    expectedResources: session.resources, actualRequest: rebuilt.renderRequest, actualResources: rebuilt.resources }) };
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
export const createGallerySessionPublication = (owners) => {
  const admit = (session) => {
    const version = Number(session?.version);
    if (version === CURRENT_VERSION) return validateCurrent(session);
    if (!HISTORICAL_VERSIONS.has(version)) throw new Error(`Gallery publication supports current version 40 or historical versions 31-33/39; received ${String(session?.version)}.`);
    return validateCurrent(owners.promoteSession(session));
  };
  const rebuild = (session) => rebuildIntent(session, owners);
  const prepare = async (source) => {
    const admitted = admit(source), { config, rebuilt, equivalence } = await rebuild(admitted);
    const session = { ...admitted, config, renderRequest: rebuilt.renderRequest, resources: rebuilt.resources, webFiles: rebuilt.webFiles };
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
