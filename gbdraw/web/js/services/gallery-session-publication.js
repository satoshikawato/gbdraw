// @ts-check
import {
  createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat,
  readStoredSimilarityAlignmentResetReceipt, validateCurrentWriterActiveConfig
} from './session-active-config-contract.js';
import {
  LEGACY_LAYOUT_PREFERENCE_FIELDS, normalizeLayoutPreferences, resolveActiveLayoutPreference, restoredLayoutPreferences
} from './layout-preferences.js';
import { validateModeSliceFields } from './mode-scoped-migration.js';
import {
  circularSourcesOfModeScoped, drawingViewOf, sessionDocumentOfViews, sessionLegendEditsByRowKey
} from './session-drawings-upgrade.js';
import {
  CANONICAL_REQUEST_SCHEMA as CURRENT_REQUEST_SCHEMA,
  FLAT_DRAFT_TOP_LEVEL_FIELDS,
  GALLERY_SESSION_DRAWN_DRAWINGS,
  MODE_SCOPED_SESSION_VERSION,
  PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  PROTEIN_LOSAT_CACHE_SCHEMA,
  SESSION_DRAWINGS_LAYOUT,
  SESSION_VERSION as CURRENT_VERSION
} from '../contracts/schema.generated.js';
import { adoptCurrentSessionResources } from './session-resource-backing.js';
/** @import { FileLike } from './session-resource-backing.js' */
/** @import { CardRecordFacts } from './linear-card-windows.js' */
import { storesDrawnFrameProteinResults } from './main-session-comparison-frame.js';
import { discoverInputRecords, readCircularSourceRecords } from './record-discovery.js';
import { circularSourceInputs, inputSourceFiles } from './input-types.js';
import { cardFileRecordFacts } from './linear-card-windows.js';
import { diagnosticError } from '../utils/error-normalization.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS, classifyImportedComparisonIntent, inheritSavedComparison, inheritUnitPlans
} from './imported-comparison-intent.js';
import { canonicalAnnotationParams, normalizeAnnotationSlotParams, validateCustomTrackPlan } from './track-slot-validation.js';
const CACHE_LIMIT_BYTES = 64 * 1024 * 1024;
const ARTIFACT_FIELDS = ['results', 'editorState', 'orthogroupState', 'runMetadata', 'losatCache', 'losatDerivedCache', 'proteinIdentityManifest'];
const noComparisonPlan = () => ({ mode: 'none', defaultSource: 'losat', edges: [] });
/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/**
 * @template T
 * @param {T} value
 * @returns {T}
 */
const clone = (value) => value === undefined ? value : JSON.parse(JSON.stringify(value));
/**
 * @param {unknown} value
 * @param {string} key
 */
const has = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
/** @import { FeatureCatalog } from './feature-catalog.js' */
/** @import { Session51 } from './session-drawings-upgrade.js' */
/** @import { LinearUnitPlanner } from './imported-comparison-intent.js' */
/** @import { ActiveWebConfig } from './session-active-config-contract.js' */
/** @import { LayoutPreferences, StoredLayoutPreferences } from './layout-preferences.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */

/** @typedef {Partial<ActiveWebConfig> & Record<string, unknown>} PublicationConfig A Session 46 slice may omit a field; it reads that mode's default (plan 4.1). */
/** @typedef {PublicationConfig & Pick<ActiveWebConfig, 'form' | 'adv'>} PublishedConfig The configuration the publication rebuilds: it always has `form` and `adv`. */

/** @typedef {{ schema?: unknown, kind?: unknown, program?: unknown, idEncoding?: unknown, queryProteinSetHash?: unknown, subjectProteinSetHash?: unknown }} PublicationLosatEntry */
/** @typedef {{ config?: PublicationConfig, features?: Record<string, unknown>, ui?: { layoutPreferences?: unknown } & Record<string, unknown> }} PublicationModeSlice */
/**
 * What `owners.projectRequest` returns for a Session's request, as the publication reads it.
 * @typedef {object} PublicationProjection
 * @property {'circular' | 'linear'} mode
 * @property {PublicationConfig} [config]
 * @property {StoredLayoutPreferences | null} [layoutPreferences]
 * @property {{ linearSeqs: PublicationLinearCard[], c_depth?: unknown,
 *   circularSources?: PublicationCircularSource[], [field: string]: unknown }} files The input files; a Circular
 *   projection names its File cards (`circularSources`).
 */
/**
 * A Linear File card of a projection: the resource of each file it holds, the record it names, and what the rebuild
 * sets before it builds the request (the committed record's `cardinality`, and `cardRecord`, its file's records).
 * @typedef {{ gb: FileLike | null, gff: FileLike | null, fasta: FileLike | null, region_record_id: string,
 *   cardinality?: unknown, depth?: unknown, cardRecord?: CardRecordFacts }} PublicationLinearCard
 */
/**
 * A Circular File card of a projection: its uid and the resource of each file it holds.
 * @typedef {{ uid: string, gb: FileLike | null, gff: FileLike | null, fasta: FileLike | null }} PublicationCircularSource
 */
/**
 * A Session document at publication. Admission checks it: the current writer format (`modes`), an
 * older Session's flat Web draft (`config`, `features`), or a CLI-written Session (neither).
 * @typedef {PublicationSessionFields & Record<string, unknown>} PublicationSession
 */
/**
 * A Session's request as publication reads it. Before admission it is unvalidated and of any schema that
 * admission promotes (5 to 9); a current request (`CanonicalRenderRequest`) is one of these.
 * @typedef {object} PublicationRequest
 * @property {unknown} [schema]
 * @property {string} mode
 * @property {{ config?: unknown, colors?: { defaultColors?: unknown, defaultColorsFile?: unknown } }} [diagramOptions]
 * @property {{ kind?: unknown }[]} [comparisons]
 * @property {unknown} [output]
 */
/**
 * @typedef {object} PublicationSessionFields
 * @property {string} [format]
 * @property {number} [version]
 * @property {PublicationRequest} renderRequest
 * @property {Record<string, { kind?: unknown, encoding?: unknown, data?: unknown }>} [resources]
 * @property {Record<string, unknown>} [webFiles]
 * @property {unknown} [files]
 * @property {PublicationConfig} [config]
 * @property {Record<string, unknown>} [features]
 * @property {Record<string, PublicationModeSlice>} [modes]
 * @property {unknown} [ui]
 * @property {{ generatedBy?: unknown, fileBindings?: unknown }} [cliInvocation]
 * @property {{ featureCatalog?: FeatureCatalog | null, alignmentResetReceipt?: unknown }} [editorState]
 * @property {Record<string, unknown> | null} [orthogroupState]
 * @property {{ entries?: PublicationLosatEntry[] }} [losatCache]
 * @property {{ entries?: unknown[] }} [losatDerivedCache]
 * @property {{ schema?: unknown }} [proteinIdentityManifest]
 * @property {unknown} [otherModeResult]
 */

/**
 * The drawing a published Session holds: its source's, or its mode's for a Session 27-50 source.
 * @typedef {{ id: string, name: string }} PublishedDrawing
 */
/**
 * @param {PublicationSession | Session51} document
 * @returns {document is Session51} The document lists its drawings (Session 51).
 */
const holdsDrawings = (document) => Array.isArray(document.drawings);
/**
 * The drawing a Gallery Session shows, as the one-drawing view the publication rebuilds
 * (`drawingViewOf`), and its ID and name. A Gallery Session shows one drawing: a Session 51's only
 * drawing with a committed render. A Web Save also holds an undrawn drawing of the other mode; the
 * publication leaves it out, as it left out the other mode's slice of a Session 46-50. A Session
 * 27-50 is its own view.
 * @param {PublicationSession | Session51} document A Session file as read, before admission.
 * @returns {{ view: PublicationSession, drawing: PublishedDrawing | null }}
 */
export const galleryPublicationView = (document) => {
  if (!holdsDrawings(document)) return { view: document, drawing: null };
  const drawn = document.drawings.filter((drawing) => drawing.renderRequest !== null);
  if (drawn.length !== GALLERY_SESSION_DRAWN_DRAWINGS) {
    // The rule is Python's (`gbdraw/session_contract.py`); the refresh tool checks it too.
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  }
  const [shown] = drawn;
  // The view holds the drawing's request (not null, tested above) at the top level.
  const view = /** @type {PublicationSession} */ (/** @type {unknown} */ (drawingViewOf(document, shown)));
  return { view, drawing: { id: shown.id, name: shown.name } };
};
/**
 * The Session 51 a publication writes of its one-drawing view (`sessionDocumentOfViews`): one drawing,
 * of the source's ID and name, or of its mode's (`drawing` null, a Session 27-50 source).
 * @param {PublicationSession} view
 * @param {PublishedDrawing | null} drawing
 */
export const galleryPublicationDocument = (view, drawing) => {
  const mode = /** @type {'circular' | 'linear'} */ (view.renderRequest.mode);
  const { id, name } = drawing ?? { id: mode, name: SESSION_DRAWINGS_LAYOUT.drawingNames[mode] };
  return sessionDocumentOfViews([{ id, name, view }], id);
};

/** @param {PublicationSession['orthogroupState']} value */
const currentOrthogroupState = (value) => {
  if (!isObject(value)) return value;
  const current = clone(value);
  delete current.selectedOrthogroupAlignmentFeature;
  return current;
};
/** @param {PublicationSession} session */
const regenerableProteinCache = (session) => session.renderRequest?.comparisons?.some(
  ({ kind }) => kind === 'generatedProteinComparison') && session.proteinIdentityManifest?.schema === PROTEIN_IDENTITY_MANIFEST_SCHEMA
  && session.losatCache?.entries?.some((entry) => entry?.schema === PROTEIN_LOSAT_CACHE_SCHEMA && entry?.kind === 'raw-losat' && entry?.program === 'blastp'
    && entry?.idEncoding === 'runtime-handle-v1'
    && typeof entry?.queryProteinSetHash === 'string' && typeof entry?.subjectProteinSetHash === 'string');
/**
 * @param {PublicationSession} session
 * @param {{ limitBytes?: number }} [options]
 * @returns {PublicationSession}
 */
export const applyDerivedCachePublicationPolicy = (session, { limitBytes = CACHE_LIMIT_BYTES } = {}) => {
  const entries = session.losatDerivedCache?.entries;
  if (!entries?.length || !regenerableProteinCache(session) || new TextEncoder().encode(JSON.stringify(entries)).byteLength <= limitBytes) return session;
  return { ...session, losatDerivedCache: { ...session.losatDerivedCache, entries: [] } };
};
// The draft of the Session's mode: its slice (Session 46) or the flat draft of
// an older Session that admission promoted. Publication writes the slice of the
// Session's mode only; the other mode takes its defaults (plan 4.3).
/** @param {PublicationSession} session */
const modeSlice = (session) => {
  const slice = isObject(session?.modes) ? session.modes[session.renderRequest?.mode] : null;
  return isObject(slice) ? slice : null;
};
// The Session holds the Web draft of its mode in that mode's slice, as Session
// Load reads it (`hasStoredDraft`, services/config.js): a CLI drawing of a
// Session 51 has a slice without a draft.
/** @param {PublicationSession} session */
const holdsSliceDraft = (session) => isObject(modeSlice(session)?.config);
// A CLI-written Session has no Web configuration; publication derives it from
// the request, as Session Load does, so a Gallery Session can be built from its
// declared command.
/** @param {PublicationSession} session */
const isCliWritten = (session) => !has(session, 'config') && !holdsSliceDraft(session)
  && session?.cliInvocation?.generatedBy === 'gbdraw';
/** @param {unknown} output */
const outputEntries = (output) => {
  const entries = Array.isArray(output) ? output : [output];
  return entries.every(isObject) ? entries : null;
};
/** @param {PublicationSession} session */
const draftConfig = (session) => (holdsSliceDraft(session) ? modeSlice(session)?.config : session?.config);
/** @param {PublicationSession} session */
const draftFeatures = (session) => (holdsSliceDraft(session) ? modeSlice(session)?.features : session?.features);
/** @param {PublicationSession} session */
const validateEnvelope = (session) => {
  if (!isObject(session) || session.format !== 'gbdraw-session') throw new Error('Gallery publication requires a gbdraw-session document.'); if (Number(session.version) !== CURRENT_VERSION) throw new Error(`Gallery publication requires session version ${CURRENT_VERSION}.`);
  if (!isObject(session.renderRequest) || Number(session.renderRequest.schema) !== CURRENT_REQUEST_SCHEMA) throw new Error(`Gallery publication requires canonical renderRequest schema ${CURRENT_REQUEST_SCHEMA}.`);
  // A Gallery Session shows one diagram mode (E1: `otherModeResult` holds a second).
  if (has(session, 'otherModeResult')) throw new Error('Gallery publication requires a Session with one diagram mode.');
  return session;
};
// A Session 46 keeps the Web draft in `modes` only (plan 4.1).
/** @param {PublicationSession} session */
const validateCurrentFields = (session) => {
  // The top-level homes of the flat draft that Session 46 moved into `modes`.
  if (FLAT_DRAFT_TOP_LEVEL_FIELDS.some((field) => has(session, field))) throw new Error('Gallery publication requires the Web draft in modes (Session 46).');
  if (has(session, 'modes')) {
    const modes = session.modes;
    if (!isObject(modes) || Object.keys(modes).some((mode) => mode !== 'circular' && mode !== 'linear'))
      throw new Error('Gallery publication requires modes to hold a circular and a linear slice only.');
    Object.values(modes).forEach((slice) => validateModeSliceFields(slice));
  }
  return session;
};
/** @param {PublicationSession} session */
const validateCurrent = (session) => {
  validateEnvelope(session); validateCurrentWriterActiveConfig({ mode: session.renderRequest.mode, storedConfig: draftConfig(session) });
  return session;
};
/**
 * @param {PublicationSession} session
 * @param {PublicationProjection} projection
 * @returns {PublishedConfig}
 */
const publicationConfig = (session, projection) => {
  const projected = clone(projection.config || {}), stored = draftConfig(session) || {}, losat = createDefaultLosat();
  /** @type {PublishedConfig} */
  const config = { ...projected,
    form: { ...createDefaultForm(), ...projected.form, ...stored.form },
    adv: { ...createDefaultAdv(projection.mode), ...projected.adv, ...stored.adv },
    losat: { ...losat, ...projected.losat, ...stored.losat,
      blastn: { ...losat.blastn, ...projected.losat?.blastn, ...stored.losat?.blastn },
      blastp: { ...losat.blastp, ...projected.losat?.blastp, ...stored.losat?.blastp } },
    circularConservation: { ...createDefaultCircularConservation(), ...projected.circularConservation, ...clone(stored.circularConservation) },
    linearComparisonPlan: clone(stored.linearComparisonPlan || projected.linearComparisonPlan || noComparisonPlan())
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
/**
 * @param {PublicationSession} session
 * @param {PublishedConfig} config
 * @param {PublicationProjection} projection
 */
const publicationLayout = (session, config, projection) => {
  const mode = projection.mode;
  const slot = modeSlice(session)?.ui?.layoutPreferences;
  const ui = isObject(session.ui) ? session.ui : {};
  const layout = holdsSliceDraft(session)
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
// that mode's defaults. The rebuilt draft is split by the registry through
// Python's split (`owners.splitDraft`); a Session 46 keeps the rest of its
// slice. The layout has one home, so the older `ui` layout fields go.
/**
 * @param {PublicationSession} session
 * @param {PublishedConfig} config
 * @param {LayoutPreferences} layout
 * @param {PublicationProjection} projection
 * @param {Record<string, unknown>} features
 * @param {PublicationDraftSplit} splitDraft
 */
const publishedSession = async (session, config, layout, projection, features, splitDraft) => {
  const mode = projection.mode;
  const flat = !holdsSliceDraft(session);
  // Display drafts that carry no scope (those projected from a CLI request and
  // those of a Session 46 slice) belong to the committed mode, so the split
  // keeps them there (OV-217).
  const scoped = Array.isArray(config.recordDisplayDrafts) ? { ...config, recordDisplayDrafts: config.recordDisplayDrafts
    .map((entry) => (isObject(entry) && !has(entry, 'scope') ? { ...entry, scope: mode } : entry)) } : config;
  // A flat draft takes the per-feature edits the request was rebuilt from: a
  // CLI Session holds its visibility and label tables only in the request, and
  // Load reads a stored draft's edits from the draft.
  // The split keeps Show Depth only for a drawing with a Depth file, which it counts in `webFiles.bindings`
  // (`session_depth_source_widths`): the files the projection bound, as Load binds them. A slice draft holds no
  // bindings of its own, so without them a Web-saved Session published with Depth hidden. The published
  // `webFiles` is the rebuilt request's (`prepare`).
  const webFiles = { bindings: { c_depth: projection.files.c_depth ?? null,
    linearSeqs: projection.files.linearSeqs.map(({ depth }) => ({ depth: depth ?? null })) } };
  const draft = flat ? { ...session, features, config: scoped, ui: { ...(isObject(session.ui) ? session.ui : {}), mode }, webFiles }
    : { config: scoped, ui: { mode }, webFiles };
  // Override colors merge into the draft palette's colors, which the split reads from the draft (OV-3257).
  // Only the committed mode's slice is published, so no mode profile feeds the split.
  const split = await splitDraft(draft, { committedMode: mode });
  const slice = flat ? split.modes[mode] : { ...modeSlice(session), config: split.modes[mode].config };
  const published = flat ? split : { ...session };
  const ui = Object.fromEntries(Object.entries(isObject(published.ui) ? published.ui : {})
    .filter(([field]) => !LEGACY_LAYOUT_PREFERENCE_FIELDS.includes(field)));
  return { ...published, ui, modes: { [mode]: { ...slice, ui: { ...(slice.ui || {}), layoutPreferences: clone(layout[mode]) } } } };
};
/**
 * The Session's request read by the one reader with the Session's context
 * (OV-771: a schema-5 record keeps the cardinality the reader's promotion gives it).
 * @param {PublicationSession} session
 * @param {GallerySessionPublicationOwners['readRequest']} readRequest
 */
const publicationCanonicalRequest = (session, readRequest) => readRequest(session.renderRequest, {
  featureCatalog: session.editorState?.featureCatalog || null,
  legacyOrthogroupState: session.orthogroupState || null,
  resourceIds: new Set(isObject(session.resources) ? Object.keys(session.resources) : [])
});
/**
 * A Session whose request the file stored at an earlier schema, read as Load reads it: the request promoted by the
 * one request reader, and its alignment reset receipt read against the promoted request by Load's reader
 * (`readStoredSimilarityAlignmentResetReceipt`), which rebinds a receipt bound to the request as stored (OV-2772).
 * A receipt that binds neither is refused: publication repairs nothing.
 * @param {PublicationSession} session
 * @param {GallerySessionPublicationOwners['readRequest']} readRequest
 * @returns {Promise<PublicationSession>}
 */
const withPromotedRequest = async (session, readRequest) => {
  const renderRequest = publicationCanonicalRequest(session, readRequest);
  const editorState = session.editorState ?? {};
  const receipt = editorState.alignmentResetReceipt;
  if (receipt === null || receipt === undefined) return { ...session, renderRequest };
  const { value, stale } = await readStoredSimilarityAlignmentResetReceipt(receipt,
    { renderRequest, resources: session.resources }, session.renderRequest);
  if (stale) throw new Error('Gallery publication cannot read the Similarity alignment reset receipt: it binds neither '
    + 'the stored nor the promoted request.');
  return { ...session, renderRequest, editorState: { ...editorState, alignmentResetReceipt: value } };
};
/**
 * Save writes each Linear File card's own state into `webFiles.bindings.linearSeqs`, among it a record's typed
 * name apart from its inferred one. A published Session keeps those rows, so Load and Reset Settings read the
 * names as after that Save (OV-1701). The imported table files stay unpublished (OV-367): only the card rows are
 * kept, and Session validation rejects a row that binds a file the Session does not hold.
 * A CLI-written Session's rows are only initial values that its request overrides, so it publishes none.
 * @param {PublicationSession} session
 * @returns {{ schema: unknown, linearSeqs: unknown[] } | null}
 */
const publishedCardBindings = (session) => {
  const bindings = !isCliWritten(session) && isObject(session.webFiles) ? session.webFiles.bindings : null;
  if (!isObject(bindings) || !Array.isArray(bindings.linearSeqs)) return null;
  return { schema: bindings.schema, linearSeqs: bindings.linearSeqs };
};
/**
 * @param {PublicationSession} session
 * @param {GallerySessionPublicationOwners} owners
 */
const rebuildIntent = async (session, owners) => {
  const renderRequest = publicationCanonicalRequest(session, owners.readRequest);
  // A read-only comparison (one the Web controls cannot represent) is rebuilt
  // as Generate's Inherit builds it (D-04, OV-268): Load projects the request
  // without it, the request draws no comparison of its own, and the committed
  // comparisons are inherited.
  const classification = classifyImportedComparisonIntent({ renderRequest, resources: session.resources });
  const inherit = classification.disposition === IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY;
  const projection = owners.projectRequest({ renderRequest: inherit ? { ...renderRequest, comparisons: [] } : renderRequest,
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
  const { state, drawing, features } = owners.buildRequestState({ session: { ...session, features: draftFeatures(session) }, projection, config });
  // The rebuild reads the records of each Circular card, as Generate does before it builds its request (D-13).
  if (projection.mode === 'circular') {
    if (!filesData.circularSources) throw new Error('Gallery publication requires the File cards of a Circular projection.');
    const reads = await readCircularSourceRecords(circularSourceInputs(drawing.inputType.value, filesData.circularSources), (file) => file);
    const failed = reads.find(({ error }) => error !== null);
    if (failed) throw failed.error;
    drawing.circularRecordList.value = reads.flatMap(({ records }) => records);
  }
  // It reads the records of each Linear card too, as Generate does (`catalogCardRecord`): the request names a card's
  // record by the records its file holds (`cliRecordSelector`).
  if (projection.mode === 'linear') {
    await Promise.all(filesData.linearSeqs.map(async (sequence) => {
      const { spec, primaryFile, pairedFile } = inputSourceFiles(drawing.inputType.value, sequence);
      sequence.cardRecord = cardFileRecordFacts(await discoverInputRecords(spec, primaryFile, pairedFile), sequence.region_record_id);
    }));
  }
  const plan = projection.mode === 'linear' ? owners.resolveComparisonPlan({
    plan: inherit ? noComparisonPlan() : drawing.linearComparisonPlan, sequences: filesData.linearSeqs,
    layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
    losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat?.blastp?.mode }) : null;
  const rebuilt = owners.buildRequest({ state, drawing, filesData: inherit ? { ...filesData, linearCanonicalComparisons: [] } : filesData,
    comparisonPlanSnapshot: plan });
  if (inherit) {
    const { source } = classification;
    const units = await inheritUnitPlans({ candidate: rebuilt, source, runHelper: owners.planLinearUnits });
    inheritSavedComparison({ candidate: rebuilt, source, units });
  }
  // A batch request holds one output per record (OV-267); each keeps its committed policy.
  const committedOutputs = outputEntries(session.renderRequest.output), rebuiltOutputs = outputEntries(rebuilt.renderRequest.output);
  if (!committedOutputs || !rebuiltOutputs || committedOutputs.length !== rebuiltOutputs.length) throw new Error('Gallery publication cannot preserve committed output metadata policy.');
  rebuiltOutputs.forEach((output, index) => { output.interactiveMetadataPolicy = committedOutputs[index].interactiveMetadataPolicy; });
  // A CLI-written request carries the resolved configuration; publication
  // writes the Web's own configOverrides, so Session Load needs no Worker. The
  // refresh tool checks that the replayed figure equals the declared figure.
  const cliConfig = isCliWritten(session) && isObject(session.renderRequest.diagramOptions?.config);
  if (isObject(session.renderRequest.diagramOptions?.config) && !cliConfig) {
    rebuilt.renderRequest.diagramOptions.config = clone(session.renderRequest.diagramOptions.config); delete rebuilt.renderRequest.diagramOptions.configOverrides;
  }
  // A CLI-written request compares without its resolved configuration;
  // defaults it leaves unstated compare as absent in the equivalence (OV-398),
  // also a Depth series style the command did not state (OV-2819).
  // Its track slots compare as Load reads them (OV-1299).
  const unstatedDepth = (renderRequest.diagramOptions.depthTracks || []).map((track) => new Set(
    Object.entries(track).filter(([, value]) => value === null).map(([key]) => key)));
  /** @param {CanonicalRenderRequest} request */
  const comparable = (request) => {
    if (!cliConfig) return request;
    const diagramOptions = Object.fromEntries(Object.entries(request.diagramOptions)
      .filter(([key]) => key !== 'config' && key !== 'configOverrides'));
    const depthTracks = Array.isArray(diagramOptions.depthTracks)
      ? diagramOptions.depthTracks.map((track, index) => (isObject(track)
        ? Object.fromEntries(Object.entries(track).filter(([key]) => !unstatedDepth[index]?.has(key)))
        : track))
      : diagramOptions.depthTracks;
    return { ...request, diagramOptions: { ...diagramOptions, ...(depthTracks === undefined ? {} : { depthTracks }),
      tracks: tracksAsLoaded(diagramOptions.tracks, projection.mode, drawing.form.track_type) } };
  };
  return { config, features, layout, projection, rebuilt, equivalence: await owners.assertRequestsEquivalent({ expectedRequest: comparable(renderRequest),
    expectedResources: session.resources, actualRequest: comparable(rebuilt.renderRequest), actualResources: rebuilt.resources }) };
};
/**
 * A CLI request's track slots as Web Load reads them (OV-1299): a slot spec's
 * annotation params are text ("0", "false"), which Load reads as the editor
 * writes them, and an Axis index the command did not state is the one the
 * Web infers from the slots.
 * @param {unknown} tracks
 * @param {string} mode
 * @param {unknown} trackType
 */
const tracksAsLoaded = (tracks, mode, trackType) => {
  const slotsKey = `${mode}TrackSlots`, axisKey = `${mode}TrackAxisIndex`;
  const slots = isObject(tracks) ? tracks[slotsKey] : null;
  if (!isObject(tracks) || !Array.isArray(slots) || slots.length === 0) return tracks;
  const read = slots.map((slot, index) => {
    if (!isObject(slot) || slot.renderer !== 'annotations' || !isObject(slot.params)) return slot;
    const params = { ...slot.params };
    normalizeAnnotationSlotParams(params, `${slotsKey}[${index}].params`);
    return { ...slot, params: canonicalAnnotationParams(params) };
  });
  const axis = tracks[axisKey] ?? validateCustomTrackPlan({ mode, slots: read, axisIndex: null, trackType }).emittedAxisIndex;
  return { ...tracks, [slotsKey]: read, [axisKey]: axis };
};
/**
 * @param {PublicationSession} prepared
 * @param {PublicationSession} replayed
 */
const mergeReplayResources = (prepared, replayed) => {
  for (const [id, expected] of Object.entries(prepared.resources || {})) {
    const actual = replayed.resources?.[id], fields = actual && /** @type {('kind' | 'encoding' | 'data')[]} */ (['kind', 'encoding', 'data']).filter((field) => expected[field] !== actual[field]);
    if (fields?.length) throw new Error(`Gallery replay resource collision at resources.${id}: ${fields.join(', ')}.`);
  }
  const resources = { ...(prepared.resources || {}) }, /** @type {Set<string>} */ referenced = new Set();
  /** @param {unknown} value */
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
 * What the draft split reads besides the draft.
 * @typedef {object} PublicationDraftSplitContext
 * @property {'circular' | 'linear'} committedMode The saved Result's mode. The draft's override colors merge into its
 *   palette's colors, which Python reads from the draft (`mode_split_palette_colors`).
 */
/**
 * A Web draft with a `config` split into Session 46 mode slices: the draft's other fields, its app-level settings
 * moved to the top level, and a slice per mode (empty when the mode takes no value).
 * @typedef {Record<string, unknown> & { ui?: Record<string, unknown>, modes: Record<'circular' | 'linear', PublicationModeSlice> }} PublicationSplitDraft
 */
/**
 * The one split of a flat Web draft into mode slices: Python's `split_draft_into_modes`
 * (`gbdraw/session_io.py`), which the publication bridge (`tools/publish_gallery_session.mjs`) runs.
 * @typedef {(draft: Record<string, unknown>, context: PublicationDraftSplitContext) => Promise<PublicationSplitDraft>} PublicationDraftSplit
 */
/**
 * The owners a publication calls: the Session and request owners
 * (`services/session-request.js`), the draft split, the Linear comparison
 * planner, and Python's plan of the drawn units (`planLinearUnits`), which
 * Inherit carries a saved comparison by. Each takes and returns unvalidated documents.
 * @typedef {object} GallerySessionPublicationOwners
 * @property {(input: Record<string, any>) => Promise<Record<string, any>>} assertRequestsEquivalent
 * @property {(input: Record<string, any>) => { renderRequest: CanonicalRenderRequest, resources: NonNullable<PublicationSession['resources']>, webFiles: Record<string, unknown> }} buildRequest
 * @property {(input: Record<string, any>) => { state: Record<string, any>, drawing: Record<string, any>, features: Record<string, unknown> }} buildRequestState
 *   The request inputs of a Session: project inputs and artifacts (`state`), settings and edits (`drawing`),
 *   and the drawing's per-feature edits as the draft stores them (`features`).
 * @property {(request: unknown, context: { featureCatalog: unknown, legacyOrthogroupState: unknown,
 *   resourceIds: ReadonlySet<string> }) => CanonicalRenderRequest} readRequest The one reader of a stored
 *   request (`readPublishedRenderRequest`), refusing any repair.
 * @property {(input: Record<string, any>) => PublicationProjection} projectRequest
 * @property {(input: Record<string, any>) => any} resolveComparisonPlan
 * @property {LinearUnitPlanner} planLinearUnits
 * @property {PublicationDraftSplit} splitDraft
 */

/**
 * @param {GallerySessionPublicationOwners} owners
 */
export const createGallerySessionPublication = (owners) => {
  /** @param {PublicationSession} session */
  const admit = async (session) => {
    const version = Number(session?.version);
    // Only Python converts the protein results these Sessions stored in the drawn frame (OV-1145).
    if (storesDrawnFrameProteinResults(session)) {
      throw new Error(`Gallery publication cannot convert the protein results that Session ${version} stored for a `
        + 'reverse-complemented record; load it in the app and save it again, or refresh it with '
        + 'tools/refresh_gallery_sessions.py.');
    }
    if (version >= MODE_SCOPED_SESSION_VERSION && version <= CURRENT_VERSION) {
      // A Session 46 differs only in the keys of its Legend edits, 47 in its protein result frame; a Session 46-50
      // binds its Circular input in `c_gb`, `c_gff` and `c_fasta`, which become File cards as at Load. A request of an
      // earlier schema (saved before windows) is read as Load reads it, promoted by the one request reader
      // (`readRequest`); the publication writes the current schema only.
      const legacy = sessionLegendEditsByRowKey(session, version).session;
      // The migration changes only the bindings and the Circular slice; `validateCurrentFields` checks the result.
      const stored = version === CURRENT_VERSION ? session : /** @type {PublicationSession} */ (
        { ...legacy, ...circularSourcesOfModeScoped(legacy), version: CURRENT_VERSION });
      const current = isObject(stored.renderRequest) && Number(stored.renderRequest.schema) !== CURRENT_REQUEST_SCHEMA
        ? await withPromotedRequest(stored, owners.readRequest) : stored;
      validateCurrentFields(current);
      return isCliWritten(current) ? validateEnvelope(current) : validateCurrent(current);
    }
    throw new Error(`Gallery publication supports Sessions ${MODE_SCOPED_SESSION_VERSION}-${CURRENT_VERSION}; Load a Session 27-44 in the app and save it (Python converts it at Load); received ${String(session?.version)}.`);
  };
  /** @param {PublicationSession} session */
  const rebuild = (session) => rebuildIntent(session, owners);
  /** @param {PublicationSession} source */
  const prepare = async (source) => {
    const admitted = await admit(source), { config, features, layout, projection, rebuilt, equivalence } = await rebuild(admitted);
    const cards = publishedCardBindings(admitted);
    const session = {
      ...(await publishedSession(admitted, config, layout, projection, features, owners.splitDraft)),
      renderRequest: rebuilt.renderRequest,
      resources: rebuilt.resources,
      webFiles: cards ? { ...rebuilt.webFiles, bindings: cards } : rebuilt.webFiles,
      orthogroupState: currentOrthogroupState(admitted.orthogroupState)
    };
    return { session: validateCurrent(validateCurrentFields(session)), equivalence };
  };
  /** @param {{ prepared: PublicationSession, replayed: PublicationSession }} input */
  const finalize = async ({ prepared, replayed }) => {
    validateCurrent(validateCurrentFields(prepared)); validateCurrent(validateCurrentFields(replayed));
    const resources = mergeReplayResources(prepared, replayed);
    await owners.assertRequestsEquivalent({ expectedRequest: prepared.renderRequest, expectedResources: prepared.resources, actualRequest: replayed.renderRequest,
      actualResources: replayed.resources });
    /** @type {PublicationSession} */
    let session = { ...prepared, resources };
    for (const field of ARTIFACT_FIELDS) if (has(replayed, field)) session[field] = replayed[field];
    // The replay decides the legacy protein artifacts, also when it writes none (OV-1144).
    if (has(replayed, 'legacyArtifacts')) session.legacyArtifacts = replayed.legacyArtifacts;
    else delete session.legacyArtifacts;
    session = applyDerivedCachePublicationPolicy(session);
    return { session, equivalence: (await rebuild(session)).equivalence };
  };
  /** @param {PublicationSession} session */
  const validate = (session) => (validateCurrent(validateCurrentFields(session)), rebuild(session));
  return Object.freeze({ admitGallerySession: admit, finalizeGallerySessionPublication: finalize, prepareGallerySessionForPublication: prepare, validateGalleryPublicationReadiness: validate });
};
