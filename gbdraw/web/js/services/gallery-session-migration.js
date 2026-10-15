// @ts-check
import {
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec
} from './circular-track-slot-model.js';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames
} from './current-option-values.js';
import {
  buildCanonicalRenderRequest,
  buildCanonicalRequestState,
  projectCanonicalSessionRequest,
  readStoredRenderRequest
} from './session-request.js';
import { migrateLegacyLinearLabelVisibility } from './linear-label-visibility.js';
import { resolveActiveLayoutPreference } from './layout-preferences.js';
import {
  createLinearComparisonEdge,
  normalizeLinearComparisonPlan,
  resolveLinearComparisonPlan
} from './linear-comparisons.js';
import { textToBase64, textToBytes } from './file-content-cache.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { SESSION_VERSION } from '../contracts/schema.generated.js';
/** @import { ActiveWebConfig } from './session-active-config-contract.js' */
/** @import { StoredLayoutPreferences } from './layout-preferences.js' */
/** @typedef {{ losatCache?: { entries?: unknown } | null }} LosatCacheSession A Session as the conservation restore reads it. */
/** @typedef {{ filename?: unknown, text?: unknown }} CachedLosatEntry */
/** @import { LinearComparisonEdge } from './linear-comparisons.js' */
/** @import { LoadRepair } from '../utils/load-repairs.js' */


/**
 * @template T
 * @param {T} value
 * @returns {value is T & Record<string, unknown>}
 */
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

/** @param {unknown} value */
const cloneJson = (value) => JSON.parse(JSON.stringify(value));

/**
 * @param {readonly unknown[]} args
 * @param {readonly string[]} names
 */
const optionValues = (args, names) => {
  const aliases = new Set(names);
  const values = [];
  for (let index = 0; index < args.length; index += 1) {
    const token = String(args[index]);
    if (aliases.has(token)) {
      if (index + 1 < args.length) values.push(String(args[index + 1]));
      index += 1;
      continue;
    }
    for (const alias of aliases) {
      if (token.startsWith(`${alias}=`)) {
        values.push(token.slice(alias.length + 1));
        break;
      }
    }
  }
  return values;
};

/** @param {{ cliInvocation?: { args?: unknown } | null, config?: { cliOptions?: { rawArgs?: unknown } | null } | null }} session */
const sessionArgs = (session) => {
  if (Array.isArray(session?.cliInvocation?.args)) {
    return session.cliInvocation.args.map((value) => String(value));
  }
  if (Array.isArray(session?.config?.cliOptions?.rawArgs)) {
    return session.config.cliOptions.rawArgs.map((value) => String(value));
  }
  return [];
};

/** @param {readonly unknown[]} args */
const explicitRepeatRendering = (args) => {
  const assignments = optionValues(
    args,
    ['--feature_shape', '--feature-shape']
  );
  for (const assignment of assignments) {
    const separator = assignment.indexOf('=');
    if (separator < 0) continue;
    const featureType = assignment.slice(0, separator).trim();
    if (featureType !== 'repeat_region') continue;
    const rendering = assignment.slice(separator + 1).trim().toLowerCase();
    return rendering || null;
  }
  return null;
};

/**
 * @param {readonly unknown[]} args
 * @param {number} index
 * @param {{ label?: string | null, subtitle?: string | null, reverseComplement?: unknown, gridRow?: number | null, gridColumn?: number | null }} [current]
 */
const presentationFromArgs = (args, index, current = {}) => {
  const labels = optionValues(args, ['--record_label', '--record-label']);
  const subtitles = optionValues(args, ['--record_subtitle', '--record-subtitle']);
  const presentation = {
    label: current?.label ?? null,
    subtitle: current?.subtitle ?? null,
    reverseComplement: Boolean(current?.reverseComplement),
    gridRow: current?.gridRow ?? null,
    gridColumn: current?.gridColumn ?? null
  };
  if (index < labels.length) presentation.label = labels[index].trim() || null;
  if (index < subtitles.length) presentation.subtitle = subtitles[index].trim() || null;
  return presentation;
};

/**
 * @param {{ records?: unknown } | null | undefined} renderRequest
 * @param {readonly unknown[]} args
 */
const hydrateRecordPresentations = (renderRequest, args) => {
  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  records.forEach((record, index) => {
    record.presentation = presentationFromArgs(args, index, record.presentation);
  });
};

/**
 * @param {{ linearSeqs?: unknown } | null | undefined} filesData
 * @param {readonly unknown[]} args
 */
const hydrateLinearFilePresentations = (filesData, args) => {
  const records = Array.isArray(filesData?.linearSeqs) ? filesData.linearSeqs : [];
  const labels = optionValues(args, ['--record_label', '--record-label']);
  const subtitles = optionValues(args, ['--record_subtitle', '--record-subtitle']);
  for (const values of [labels, subtitles]) {
    if (values.length > records.length) {
      throw new Error('CLI record presentation count exceeds the canonical record count.');
    }
  }
  records.forEach((record, index) => {
    if (index < labels.length) record.definition = labels[index];
    if (index < subtitles.length) record.record_subtitle = subtitles[index];
  });
};

/**
 * @param {CachedLosatEntry | null | undefined} entry
 * @param {string} [fallbackName]
 */
const cachedLosatFile = (entry, fallbackName) => {
  const text = String(entry?.text || '');
  const bytes = textToBytes(text);
  return {
    name: String(entry?.filename || fallbackName || 'comparison.tsv'),
    type: 'text/tab-separated-values',
    size: bytes.byteLength,
    lastModified: 0,
    encoding: 'base64',
    data: textToBase64(text)
  };
};

/** @param {LosatCacheSession | null | undefined} session */
const conservationCacheEntries = (session) => (
  (Array.isArray(session?.losatCache?.entries) ? session.losatCache.entries : [])
    .filter((entry) => (
      isPlainObject(entry) &&
      entry.kind === 'raw-losat' &&
      String(entry.program || '').toLowerCase() === 'blastn' &&
      entry.flow === 'circular-conservation' &&
      entry.display !== false &&
      typeof entry.text === 'string' &&
      entry.text.length > 0
    ))
);

/** @param {{ filename?: unknown } | null | undefined} entry */
const conservationEntryLabel = (entry) => (
  String(entry?.filename || '')
    .replace(/\.circular_conservation\.losatn\.tsv$/i, '')
    .replace(/\.losatn\.tsv$/i, '')
);

/**
 * @param {LosatCacheSession} session
 * @param {Record<string, unknown>} filesData
 * @param {{ enabled?: unknown, series?: unknown } | null | undefined} circularConservation
 */
const restoreConservationFiles = (session, filesData, circularConservation) => {
  if (Array.isArray(filesData.c_conservation_blasts) && filesData.c_conservation_blasts.length) {
    return;
  }
  if (!circularConservation?.enabled) return;
  const series = Array.isArray(circularConservation.series)
    ? circularConservation.series
    : [];
  const entries = conservationCacheEntries(session);
  if (series.length === 0 && entries.length === 0) return;
  if (series.length !== entries.length) {
    throw new Error(
      `Circular conservation has ${series.length} series but ${entries.length} reusable LOSAT result(s).`
    );
  }
  const remaining = new Set(entries.map((_, index) => index));
  filesData.c_conservation_blasts = series.map((item, seriesIndex) => {
    const label = String(item?.label || '').trim();
    if (!label) {
      throw new Error(`Circular conservation series #${seriesIndex + 1} has no label.`);
    }
    const entryIndex = entries.findIndex((entry, index) => (
      remaining.has(index) && label && conservationEntryLabel(entry) === label
    ));
    if (entryIndex < 0) {
      throw new Error(`Circular conservation result is missing for series '${label}'.`);
    }
    remaining.delete(entryIndex);
    return cachedLosatFile(
      entries[entryIndex],
      `${label || `comparison-${seriesIndex + 1}`}.circular_conservation.losatn.tsv`
    );
  });
  // Preserve these derived tables as typed-render inputs. Ordinary LOSAT saves
  // must not mistake stale files from the upload mode for current results.
  filesData.c_conservation_blasts_source = 'losat-cache';
};

/**
 * @param {{ config?: Partial<ActiveWebConfig> }} session
 * @param {{ config?: Partial<ActiveWebConfig>, layoutPreferences: StoredLayoutPreferences | null | undefined, mode: string }} projection
 */
const mergedGuiConfig = (session, projection) => {
  const saved = migratePersistedGalleryConfig(
    isPlainObject(session.config) ? session.config : {}
  );
  const projected = isPlainObject(projection.config) ? projection.config : {};
  const committedLayout = resolveActiveLayoutPreference(
    projection.layoutPreferences,
    projection.mode,
    Boolean(projected.form?.multi_record_canvas)
  );
  return {
    ...saved,
    ...projected,
    form: { ...(saved.form || {}), ...(projected.form || {}), legend: committedLayout.legend },
    adv: {
      ...(saved.adv || {}),
      ...(projected.adv || {}),
      plot_title_position: committedLayout.plotTitlePosition
    },
    palette: String(projected.palette || 'default'),
    colors: cloneJson(isPlainObject(projected.colors) ? projected.colors : {}),
    annotationSets: Array.isArray(projected.annotationSets)
      ? cloneJson(projected.annotationSets)
      : cloneJson(saved.annotationSets || []),
    circularConservation: isPlainObject(projected.circularConservation)
      ? cloneJson(projected.circularConservation)
      : cloneJson(isPlainObject(saved.circularConservation) ? saved.circularConservation : {})
  };
};

/** @param {Partial<ActiveWebConfig>} config */
const migratePersistedGalleryConfig = (config) => {
  const migratedNames = migratePersistedWebStateFieldNames(config);
  const form = { ...(migratedNames.form || {}) };
  const adv = { ...(migratedNames.adv || {}) };
  form.linear_track_layout = migratePersistedLinearTrackLayout(
    form.linear_track_layout
  );
  adv.label_placement = migratePersistedLinearLabelPlacement(
    adv.label_placement
  );
  adv.multi_record_size_mode = migratePersistedCircularMultiRecordSizeMode(
    adv.multi_record_size_mode
  );
  if (Array.isArray(adv.circular_track_slots)) {
    adv.circular_track_slots = adv.circular_track_slots.map((slot, index) => (
      typeof slot === 'string'
        ? parseCircularTrackSlotSpec(
            migrateLegacyCircularTrackSlotSpec(slot),
            index,
            adv.nt || 'GC',
            form.track_type || 'tuckin'
          )
        : migrateLegacyCircularTrackSlot(slot)
    ));
  }
  return {
    ...migratedNames,
    form,
    adv: migrateLegacyLinearLabelVisibility(adv, { load: true })
  };
};

/**
 * @param {unknown} value
 * @param {'upload' | 'losat'} [fallback]
 */
const legacyComparisonSource = (value, fallback = 'losat') => {
  const source = String(value || '').trim().toLowerCase();
  if (['upload', 'files', 'file'].includes(source)) return 'upload';
  if (source === 'losat') return 'losat';
  return fallback;
};

/**
 * @param {string} prefix
 * @param {number} index
 * @param {unknown} queryUid
 * @param {unknown} subjectUid
 */
const stableMigratedComparisonId = (prefix, index, queryUid, subjectUid) => {
  /** @param {unknown} value */
  const safe = (value) => String(value || '')
    .replace(/[^A-Za-z0-9._-]+/g, '-')
    .replace(/^-+|-+$/g, '') || 'record';
  return `linear-comparison-migrated-${prefix}-${index + 1}-${safe(queryUid)}-${safe(subjectUid)}`;
};

/** @param {Partial<ActiveWebConfig> & { blastSource?: unknown }} config The legacy draft: `blastSource` is the pre-plan comparison source. */
const stripLegacyComparisonConfig = (config) => {
  delete config.blastSource;
  if (isPlainObject(config.adv)) delete config.adv.blastSource;
  if (isPlainObject(config.linearRecordLayout)) {
    delete config.linearRecordLayout.comparisons;
  }
};

/** @param {{ file?: unknown } | null | undefined} comparison */
const legacyComparisonFile = (comparison) => comparison?.file || null;

/**
 * @param {{ config?: Record<string, any>, filesData?: Record<string, any>, forceWebDraft?: boolean | null }} [draft]
 * An unvalidated active configuration and its legacy file state.
 * @returns {{ config: Record<string, any>, filesData: Record<string, any> }}
 */
export const migrateLegacyLinearComparisonDraft = ({
  config: sourceConfig,
  filesData: sourceFiles,
  forceWebDraft = null
} = {}) => {
  const config = cloneJson(isPlainObject(sourceConfig) ? sourceConfig : {});
  /** @type {Record<string, any>} */
  const filesData = {
    ...(isPlainObject(sourceFiles) ? sourceFiles : {}),
    linearSeqs: (Array.isArray(sourceFiles?.linearSeqs) ? sourceFiles.linearSeqs : [])
      .map((/** @type {Record<string, unknown>} */ sequence) => ({ ...sequence }))
  };
  const legacyRows = filesData.linearSeqs.map((/** @type {Record<string, unknown>} */ sequence) => ({
    uid: String(sequence?.uid || ''),
    blast: sequence?.blast || null,
    losatFilename: String(sequence?.losat_filename || '')
  }));
  const legacyFileComparisons = (Array.isArray(sourceFiles?.linearComparisons)
    ? sourceFiles.linearComparisons
    : []).map((/** @type {Record<string, unknown>} */ comparison) => ({ ...comparison }));
  const layout = isPlainObject(config.linearRecordLayout)
    ? config.linearRecordLayout
    : null;
  const explicitComparisons = Array.isArray(layout?.comparisons)
    ? layout.comparisons.map((/** @type {Record<string, unknown>} */ comparison) => ({ ...comparison }))
    : null;
  const hasWebDraft = forceWebDraft === null
    ? Boolean(layout) || !isPlainObject(config.cliOptions)
    : forceWebDraft === true;
  const rawGlobalSource = config.blastSource ?? config.adv?.blastSource;
  const globalSource = legacyComparisonSource(rawGlobalSource);
  const legacyNone = String(rawGlobalSource || '').trim().toLowerCase() === 'none';
  // Legacy CLI writers selected LOSATP exactly when --protein_blastp_mode was
  // not `none`; that run compared adjacent records.
  const cliProteinMode = String(optionValues(
    Array.isArray(config.cliOptions?.rawArgs) ? config.cliOptions.rawArgs : [],
    ['--protein_blastp_mode', '--protein-blastp-mode']
  ).at(-1) || 'none').trim().toLowerCase();

  stripLegacyComparisonConfig(config);
  filesData.linearSeqs = filesData.linearSeqs.map((/** @type {Record<string, unknown>} */ sequence) => {
    const { blast: _blast, losat_filename: _losatFilename, ...rest } = sequence;
    return rest;
  });

  if (isPlainObject(config.linearComparisonPlan)) {
    const fileById = new Map(
      /** @type {[string, any][]} */ (legacyFileComparisons
        .map((comparison) => [String(comparison?.id || ''), legacyComparisonFile(comparison)])
        .filter(([id, file]) => id && file))
    );
    const normalized = normalizeLinearComparisonPlan(config.linearComparisonPlan);
    config.linearComparisonPlan = {
      mode: normalized.mode,
      defaultSource: normalized.defaultSource,
      edges: normalized.edges.map(({ file: _file, ...metadata }) => metadata)
    };
    filesData.linearComparisons = normalized.edges
      .map((edge) => ({ id: edge.id, file: edge.file || fileById.get(edge.id) || null }))
      .filter((binding) => binding.file);
    return { config, filesData };
  }

  if (!hasWebDraft) {
    // A CLI-only draft gains no synthetic Web comparison (No comparison). Its
    // own protein comparison is stated as the adjacent LOSATP plan.
    delete config.linearRecordLayout;
    if (cliProteinMode !== 'none') {
      config.linearComparisonPlan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
    } else {
      delete config.linearComparisonPlan;
    }
    filesData.linearComparisons = [];
    return { config, filesData };
  }

  const uidIndex = new Map(legacyRows.map((/** @type {{ uid: string }} */ row, /** @type {number} */ index) => [row.uid, index]));
  const legacyBindingEntries = legacyFileComparisons.map((comparison, /** @type {number} */ index) => ({
    index,
    comparison,
    id: String(comparison?.id || ''),
    queryUid: String(comparison?.queryUid || ''),
    subjectUid: String(comparison?.subjectUid || ''),
    file: legacyComparisonFile(comparison)
  }));
  const fileById = new Map();
  const fileByPair = new Map();
  legacyBindingEntries.forEach((entry) => {
    const { id, queryUid, subjectUid, file } = entry;
    if (id && file && !fileById.has(id)) fileById.set(id, entry);
    if (queryUid && subjectUid && file) {
      const key = `${queryUid}\u001f${subjectUid}`;
      if (!fileByPair.has(key)) fileByPair.set(key, entry);
    }
  });
  const consumedBindingIndexes = new Set();
  /** @param {{ file?: unknown, index: number } | null | undefined} entry */
  const consumeBinding = (entry) => {
    if (!entry?.file) return null;
    consumedBindingIndexes.add(entry.index);
    return entry.file;
  };
  /** @param {number} index */
  const positionalFile = (index) => {
    const row = legacyRows[index];
    const next = legacyRows[index + 1];
    if (!row || !next) return null;
    const endpointBinding = fileByPair.get(`${row.uid}\u001f${next.uid}`);
    if (row.blast) {
      // Canonical projection mirrored adjacent uploads into both legacy shapes.
      consumeBinding(endpointBinding);
      return row.blast;
    }
    return consumeBinding(endpointBinding);
  };
  const usedPayloadGaps = new Set();
  const usedIds = new Set();
  /**
   * @param {unknown} candidate
   * @param {string} prefix
   * @param {number} index
   * @param {unknown} queryUid
   * @param {unknown} subjectUid
   */
  const uniqueId = (candidate, prefix, index, queryUid, subjectUid) => {
    const base = String(candidate || '').trim()
      || stableMigratedComparisonId(prefix, index, queryUid, subjectUid);
    let id = base;
    let suffix = 2;
    while (usedIds.has(id)) {
      id = `${base}-${suffix}`;
      suffix += 1;
    }
    usedIds.add(id);
    return id;
  };
  /** @type {LinearComparisonEdge[]} */
  const edges = [];
  /** @param {Record<string, any>} edge */
  const addEdge = ({
    id,
    queryUid,
    subjectUid,
    source,
    included,
    file,
    fileActive,
    losatFilename,
    losatFilenameActive,
    prefix,
    index
  }) => {
    if (!queryUid || !subjectUid) return;
    edges.push(createLinearComparisonEdge({
      id: uniqueId(id, prefix, index, queryUid, subjectUid),
      queryUid,
      subjectUid,
      source,
      included,
      file,
      fileActive,
      losatFilename,
      losatFilenameActive
    }));
  };

  const authoritativeExplicit = Boolean(layout?.enabled) && explicitComparisons !== null;
  let mode = legacyNone ? 'none' : 'adjacent';
  if (authoritativeExplicit) {
    mode = explicitComparisons.length > 0 ? 'selected' : 'none';
    explicitComparisons.forEach((/** @type {Record<string, unknown>} */ comparison, /** @type {number} */ index) => {
      const queryIndex = Number.isInteger(Number(comparison?.queryIndex))
        ? Number(comparison.queryIndex)
        : uidIndex.get(String(comparison?.queryUid || ''));
      const subjectIndex = Number.isInteger(Number(comparison?.subjectIndex))
        ? Number(comparison.subjectIndex)
        : uidIndex.get(String(comparison?.subjectUid || ''));
      const queryUid = String(
        comparison?.queryUid || legacyRows[queryIndex]?.uid || ''
      );
      const subjectUid = String(
        comparison?.subjectUid || legacyRows[subjectIndex]?.uid || ''
      );
      const adjacentGap = Number.isInteger(queryIndex)
        && subjectIndex === queryIndex + 1
        ? queryIndex
        : -1;
      const file = consumeBinding(fileById.get(String(comparison?.id || '')))
        || consumeBinding(fileByPair.get(`${queryUid}\u001f${subjectUid}`))
        || (adjacentGap >= 0 ? positionalFile(adjacentGap) : null);
      const losatFilename = adjacentGap >= 0
        ? legacyRows[adjacentGap]?.losatFilename || ''
        : '';
      if (adjacentGap >= 0 && (file || losatFilename)) usedPayloadGaps.add(adjacentGap);
      addEdge({
        id: comparison?.id,
        queryUid,
        subjectUid,
        source: legacyComparisonSource(comparison?.source, globalSource),
        included: true,
        file,
        fileActive: Boolean(file),
        losatFilename,
        losatFilenameActive: Boolean(losatFilename),
        prefix: 'selected',
        index
      });
    });
  }

  legacyRows.slice(0, -1).forEach((/** @type {{ uid: string, losatFilename: string }} */ row, /** @type {number} */ index) => {
    const file = positionalFile(index);
    const losatFilename = row.losatFilename;
    if ((!file && !losatFilename) || usedPayloadGaps.has(index)) return;
    const payloadActive = mode === 'adjacent';
    const fileActive = payloadActive && globalSource === 'upload' && Boolean(file);
    const losatFilenameActive = payloadActive && globalSource === 'losat'
      && Boolean(losatFilename);
    addEdge({
      queryUid: row.uid,
      subjectUid: legacyRows[index + 1].uid,
      source: globalSource,
      included: fileActive || losatFilenameActive,
      file,
      fileActive,
      losatFilename,
      losatFilenameActive,
      prefix: 'adjacent',
      index
    });
  });

  legacyBindingEntries.forEach((entry) => {
    if (!entry.file || consumedBindingIndexes.has(entry.index)) return;
    const queryIndex = Number.isInteger(Number(entry.comparison?.queryIndex))
      ? Number(entry.comparison.queryIndex)
      : uidIndex.get(entry.queryUid);
    const subjectIndex = Number.isInteger(Number(entry.comparison?.subjectIndex))
      ? Number(entry.comparison.subjectIndex)
      : uidIndex.get(entry.subjectUid);
    const queryUid = entry.queryUid || legacyRows[queryIndex]?.uid || '';
    const subjectUid = entry.subjectUid || legacyRows[subjectIndex]?.uid || '';
    addEdge({
      id: entry.id,
      queryUid,
      subjectUid,
      source: legacyComparisonSource(entry.comparison?.source, globalSource),
      included: false,
      file: entry.file,
      fileActive: false,
      losatFilename: '',
      losatFilenameActive: false,
      prefix: 'retained',
      index: entry.index
    });
  });

  const planWithFiles = normalizeLinearComparisonPlan({
    mode,
    defaultSource: globalSource,
    edges
  });
  config.linearComparisonPlan = {
    mode: planWithFiles.mode,
    defaultSource: planWithFiles.defaultSource,
    edges: planWithFiles.edges.map(({ file: _file, ...metadata }) => metadata)
  };
  filesData.linearComparisons = planWithFiles.edges
    .filter((edge) => edge.file)
    .map((edge) => ({ id: edge.id, file: edge.file }));
  return { config, filesData };
};

/** @typedef {{ kind?: unknown, resourceId?: unknown, settings?: { alignOrthogroupFeature?: unknown, [field: string]: unknown } | null }} PreservedComparison */
/**
 * @param {{ renderRequest?: { comparisons?: unknown } | null, resources?: Record<string, unknown> | null }} session
 * @param {{ renderRequest: { comparisons?: PreservedComparison[] }, resources: Record<string, unknown> }} promoted
 */
const preserveComparisonResources = (session, promoted) => {
  const comparisons = Array.isArray(session?.renderRequest?.comparisons)
    ? cloneJson(session.renderRequest.comparisons)
    : [];
  const preservingSavedComparisons = comparisons.length > 0;
  if (preservingSavedComparisons) promoted.renderRequest.comparisons = comparisons;
  for (const comparison of promoted.renderRequest.comparisons || []) {
    if (comparison?.kind === 'generatedProteinComparison' && comparison.settings) {
      delete comparison.settings.alignOrthogroupFeature;
    }
  }
  promoted.resources = {
    ...(session.resources || {}),
    ...promoted.resources
  };
  for (const comparison of promoted.renderRequest.comparisons || []) {
    const resourceId = String(comparison?.resourceId || '');
    if (resourceId && preservingSavedComparisons && session.resources?.[resourceId]) {
      promoted.resources[resourceId] = session.resources[resourceId];
    }
    if (resourceId && !promoted.resources[resourceId]) {
      throw new Error(`Canonical comparison resource is missing: ${resourceId}`);
    }
  }
};

/**
 * @param {Record<string, any>} session A validated Session of version 31-33 or 39.
 * @param {string[]} args
 * @returns {Record<string, any>} A current-request Session with its flat draft (`config`), which publication and Load split into mode slices.
 */
const promoteCliAuthoredSession = (session, args) => {
  const sourceConfig = session.renderRequest.diagramOptions?.config;
  const promoted = promoteGuiAuthoredSession(session, args, false);
  const renderRequest = promoted.renderRequest;
  const featureShapes = isPlainObject(renderRequest.diagramOptions.featureShapes)
    ? renderRequest.diagramOptions.featureShapes
    : {};
  const explicitRepeat = explicitRepeatRendering(args);
  featureShapes.repeat_region = explicitRepeat || 'underlay';
  renderRequest.diagramOptions.featureShapes = featureShapes;
  if (sourceConfig !== undefined) {
    renderRequest.diagramOptions.config = cloneJson(sourceConfig);
  }
  return promoted;
};

/**
 * @param {Record<string, any>} session A validated Session of version 31-33 or 39.
 * @param {string[]} args
 * @param {boolean} [forceWebDraft]
 * @returns {Record<string, any>} A current-request Session with its flat draft (`config`), which publication and Load split into mode slices.
 */
const promoteGuiAuthoredSession = (session, args, forceWebDraft = true) => {
  const projection = projectCanonicalSessionRequest({
    renderRequest: session.renderRequest,
    resources: session.resources,
    webFiles: session.webFiles || {},
    legacyFiles: session.files,
    storedConfig: session.config,
    fileBindings: session.cliInvocation?.fileBindings
  });
  const migratedConfig = migratePersistedGalleryConfig(
    mergedGuiConfig(session, projection)
  );
  const projectedFiles = {
    ...projection.files,
    linearSeqs: (projection.files.linearSeqs || []).map((/** @type {Record<string, unknown>} */ record) => ({ ...record })),
    linearComparisons: (projection.files.linearComparisons || []).map((/** @type {Record<string, unknown>} */ comparison) => ({
      ...comparison
    }))
  };
  const migratedDraft = migrateLegacyLinearComparisonDraft({
    config: migratedConfig,
    filesData: projectedFiles,
    forceWebDraft
  });
  const config = migratedDraft.config;
  const filesData = migratedDraft.filesData;
  hydrateLinearFilePresentations(filesData, args);
  const { state, drawing } = buildCanonicalRequestState({ session, projection, config });
  restoreConservationFiles(session, filesData, drawing.circularConservation);
  const comparisonPlanSnapshot = projection.mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan,
        sequences: filesData.linearSeqs,
        layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
        losatProgram: drawing.losatProgram.value,
        blastpMode: drawing.losat?.blastp?.mode
      })
    : null;
  const promotedCore = buildCanonicalRenderRequest({
    state,
    drawing,
    filesData,
    comparisonPlanSnapshot
  });
  const orthogroupState = isPlainObject(session.orthogroupState)
    ? cloneJson(session.orthogroupState)
    : session.orthogroupState;
  if (isPlainObject(orthogroupState)) {
    delete orthogroupState.selectedOrthogroupAlignmentFeature;
  }
  const promoted = ({
    ...session,
    format: 'gbdraw-session',
    // The flat draft goes to its mode slice when publication prepares it.
    version: SESSION_VERSION,
    config: cloneJson(config),
    renderRequest: promotedCore.renderRequest,
    resources: promotedCore.resources,
    orthogroupState,
    webFiles: {
      ...(session.webFiles || {}),
      ...promotedCore.webFiles
    }
  });
  hydrateRecordPresentations(promoted.renderRequest, args);
  preserveComparisonResources(session, promoted);
  return promoted;
};

/**
 * A historical Gallery Session whose request its reader has read (Load reads
 * it with the Session's other requests) as a current-request Session.
 * @param {Parameters<typeof promoteGuiAuthoredSession>[0]} session A Session of version 31-33 or 39 with a read request.
 * @returns {ReturnType<typeof promoteGuiAuthoredSession>} A current-request Session with its flat draft (`config`).
 */
export const promoteReadGallerySession = (session) => {
  const args = sessionArgs(session);
  const cliAuthored = isPlainObject(session?.config?.cliOptions)
    && !isPlainObject(session?.config?.linearRecordLayout);
  return cliAuthored
    ? promoteCliAuthoredSession(session, args)
    : promoteGuiAuthoredSession(session, args);
};

/**
 * A historical Gallery Session as a current-request Session: its request is
 * read by the one reader (`readStoredRenderRequest`, which lifts schemas 1
 * and 2), then the Session is promoted.
 * @param {unknown} session An unvalidated historical Gallery Session.
 * @param {{ repairs: LoadRepair[] }} options `repairs` receives what the read of its request could not carry.
 * @returns {Record<string, any>} A current-request Session with its flat draft (`config`), which publication and Load split into mode slices.
 */
export const promoteGallerySessionToCurrent = (session, { repairs }) => {
  if (!isPlainObject(session) || !isPlainObject(session.renderRequest)) {
    throw new Error('Gallery session must contain a canonical renderRequest.');
  }
  if (!isPlainObject(session.resources)) {
    throw new Error('Gallery session must contain canonical resources.');
  }
  if (![31, 32, 33, 39].includes(Number(session.version))) {
    throw new Error(
      `Gallery historical migration supports session versions 31-33 and 39; received ${String(session.version)}.`
    );
  }
  const read = readStoredRenderRequest(session.renderRequest, {
    featureCatalog: isPlainObject(session.editorState) ? session.editorState.featureCatalog ?? null : null,
    legacyOrthogroupState: session.orthogroupState ?? null,
    resourceIds: new Set(Object.keys(session.resources)),
    path: 'renderRequest'
  });
  if ('kind' in read) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  repairs.push(...read.repairs);
  return promoteReadGallerySession({ ...session, renderRequest: read.value });
};
