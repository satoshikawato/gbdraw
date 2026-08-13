import { normalizePaletteColors } from '../app/color-utils.js';
import {
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec
} from '../app/circular-track-slots.js';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames
} from '../app/current-option-values.js';
import {
  CANONICAL_REQUEST_SCHEMA,
  buildCanonicalRenderRequest,
  projectCanonicalSessionRequest
} from './session-request.js';
import {
  createLinearComparisonEdge,
  normalizeLinearComparisonPlan,
  resolveLinearComparisonPlan
} from '../app/linear-comparisons.js';
import { textToBase64, textToBytes } from './file-content-cache.js';

const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

const cloneJson = (value) => JSON.parse(JSON.stringify(value));
const ref = (value) => ({ value });

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

const sessionArgs = (session) => {
  if (Array.isArray(session?.cliInvocation?.args)) {
    return session.cliInvocation.args.map((value) => String(value));
  }
  if (Array.isArray(session?.config?.cliOptions?.rawArgs)) {
    return session.config.cliOptions.rawArgs.map((value) => String(value));
  }
  return [];
};

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

const hydrateRecordPresentations = (renderRequest, args) => {
  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  records.forEach((record, index) => {
    record.presentation = presentationFromArgs(args, index, record.presentation);
  });
};

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

const selectorToken = (record, index) => {
  const selector = record?.region?.selector || record?.selector;
  if (selector?.kind === 'recordId') return String(selector.value || '');
  if (selector?.kind === 'recordIndex') return `#${Number(selector.index) + 1}`;
  return `#${index + 1}`;
};

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

const conservationEntryLabel = (entry) => (
  String(entry?.filename || '')
    .replace(/\.circular_conservation\.losatn\.tsv$/i, '')
    .replace(/\.losatn\.tsv$/i, '')
);

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

const mergedGuiConfig = (session, projection) => {
  const saved = migratePersistedWebStateFieldNames(
    isPlainObject(session.config) ? session.config : {}
  );
  const projected = isPlainObject(projection.config) ? projection.config : {};
  return {
    ...saved,
    ...projected,
    form: { ...(saved.form || {}), ...(projected.form || {}) },
    adv: { ...(saved.adv || {}), ...(projected.adv || {}) },
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
  return { ...migratedNames, form, adv };
};

const legacyComparisonSource = (value, fallback = 'losat') => {
  const source = String(value || '').trim().toLowerCase();
  if (['upload', 'files', 'file'].includes(source)) return 'upload';
  if (source === 'losat') return 'losat';
  return fallback;
};

const stableMigratedComparisonId = (prefix, index, queryUid, subjectUid) => {
  const safe = (value) => String(value || '')
    .replace(/[^A-Za-z0-9._-]+/g, '-')
    .replace(/^-+|-+$/g, '') || 'record';
  return `linear-comparison-migrated-${prefix}-${index + 1}-${safe(queryUid)}-${safe(subjectUid)}`;
};

const stripLegacyComparisonConfig = (config) => {
  delete config.blastSource;
  if (isPlainObject(config.adv)) delete config.adv.blastSource;
  if (isPlainObject(config.linearRecordLayout)) {
    delete config.linearRecordLayout.comparisons;
  }
};

const legacyComparisonFile = (comparison) => comparison?.file || null;

export const migrateLegacyLinearComparisonDraft = ({
  config: sourceConfig,
  filesData: sourceFiles,
  forceWebDraft = null
} = {}) => {
  const config = cloneJson(isPlainObject(sourceConfig) ? sourceConfig : {});
  const filesData = {
    ...(isPlainObject(sourceFiles) ? sourceFiles : {}),
    linearSeqs: (Array.isArray(sourceFiles?.linearSeqs) ? sourceFiles.linearSeqs : [])
      .map((sequence) => ({ ...sequence }))
  };
  const legacyRows = filesData.linearSeqs.map((sequence) => ({
    uid: String(sequence?.uid || ''),
    blast: sequence?.blast || null,
    losatFilename: String(sequence?.losat_filename || '')
  }));
  const legacyFileComparisons = (Array.isArray(sourceFiles?.linearComparisons)
    ? sourceFiles.linearComparisons
    : []).map((comparison) => ({ ...comparison }));
  const layout = isPlainObject(config.linearRecordLayout)
    ? config.linearRecordLayout
    : null;
  const explicitComparisons = Array.isArray(layout?.comparisons)
    ? layout.comparisons.map((comparison) => ({ ...comparison }))
    : null;
  const hasWebDraft = forceWebDraft === null
    ? Boolean(layout) || !isPlainObject(config.cliOptions)
    : forceWebDraft === true;
  const rawGlobalSource = config.blastSource ?? config.adv?.blastSource;
  const globalSource = legacyComparisonSource(rawGlobalSource);
  const legacyNone = String(rawGlobalSource || '').trim().toLowerCase() === 'none';

  stripLegacyComparisonConfig(config);
  filesData.linearSeqs = filesData.linearSeqs.map((sequence) => {
    const { blast: _blast, losat_filename: _losatFilename, ...rest } = sequence;
    return rest;
  });

  if (isPlainObject(config.linearComparisonPlan)) {
    const fileById = new Map(
      legacyFileComparisons
        .map((comparison) => [String(comparison?.id || ''), legacyComparisonFile(comparison)])
        .filter(([id, file]) => id && file)
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
    delete config.linearRecordLayout;
    delete config.linearComparisonPlan;
    filesData.linearComparisons = [];
    return { config, filesData };
  }

  const uidIndex = new Map(legacyRows.map((row, index) => [row.uid, index]));
  const legacyBindingEntries = legacyFileComparisons.map((comparison, index) => ({
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
  const consumeBinding = (entry) => {
    if (!entry?.file) return null;
    consumedBindingIndexes.add(entry.index);
    return entry.file;
  };
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
  const edges = [];
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
    explicitComparisons.forEach((comparison, index) => {
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

  legacyRows.slice(0, -1).forEach((row, index) => {
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

const circularConservationState = (config) => {
  const conservation = cloneJson(
    isPlainObject(config.circularConservation)
      ? config.circularConservation
      : { enabled: false, reference: 'auto', labels: '', series: [] }
  );
  const series = Array.isArray(conservation.series) ? conservation.series : [];
  if (!String(conservation.labels || '').trim() && series.length > 0) {
    conservation.labels = series.map((entry) => String(entry?.label || '').trim()).join(',');
  }
  return conservation;
};

const buildStateFacade = (session, projection, config) => {
  const features = isPlainObject(session.features) ? session.features : {};
  const layout = isPlainObject(config.linearRecordLayout) ? config.linearRecordLayout : {};
  const palette = String(config.palette || 'default');
  return {
    mode: ref(projection.mode),
    cInputType: ref(session?.ui?.cInputType || projection.inputType),
    lInputType: ref(session?.ui?.lInputType || projection.inputType),
    circularRecordList: ref(
      (session.renderRequest.records || []).map((record, index) => ({
        selector: selectorToken(record, index)
      }))
    ),
    form: config.form || {},
    adv: config.adv || {},
    normalizePaletteColors,
    paletteDefinitions: ref({ [palette]: {} }),
    currentColors: ref(config.colors || {}),
    selectedPalette: ref(palette),
    manualSpecificRules: cloneJson(config.rules || []),
    featureVisibilityRules: ref(cloneJson(
      features.featureVisibilityManualRules ||
      projection.semanticFeatureState?.featureVisibilityManualRules ||
      []
    )),
    filterMode: ref(config.filterMode || 'None'),
    manualBlacklist: ref(String(config.blacklistText || '')),
    manualWhitelist: cloneJson(config.whitelist || []),
    manualPriorityRules: cloneJson(config.qualifierPriorityRules || []),
    labelTextFeatureOverrides: cloneJson(features.labelTextFeatureOverrides || {}),
    labelTextBulkOverrides: cloneJson(features.labelTextBulkOverrides || {}),
    labelTextFeatureOverrideSources: cloneJson(features.labelTextFeatureOverrideSources || {}),
    labelVisibilityOverrides: cloneJson(features.labelVisibilityOverrides || {}),
    canonicalLabelOverrideRows: ref(cloneJson(features.labelOverrideRows || [])),
    editableLabels: ref([]),
    extractedFeatures: ref(features.extractedFeatures || []),
    circularConservation: circularConservationState(config),
    losatProgram: ref(config.losatProgram || 'blastn'),
    losat: cloneJson(config.losat || { blastp: {} }),
    selectedOrthogroupAlignmentFeature: ref(
      session?.orthogroupState?.selectedOrthogroupAlignmentFeature || ''
    ),
    linearRecordLayoutEnabled: ref(Boolean(layout.enabled)),
    linearRecordGap: ref(layout.recordGap ?? 24),
    linearRecordRows: cloneJson(layout.rows || []),
    linearComparisonPlan: normalizeLinearComparisonPlan(
      config.linearComparisonPlan || { mode: 'none', defaultSource: 'losat', edges: [] }
    ),
    annotationSets: cloneJson(config.annotationSets || [])
  };
};

const preserveComparisonResources = (session, promoted) => {
  const comparisons = Array.isArray(session?.renderRequest?.comparisons)
    ? cloneJson(session.renderRequest.comparisons)
    : [];
  const preservingSavedComparisons = comparisons.length > 0;
  if (preservingSavedComparisons) promoted.renderRequest.comparisons = comparisons;
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

const promoteGuiAuthoredSession = (session, args, forceWebDraft = true) => {
  const projection = projectCanonicalSessionRequest({
    renderRequest: session.renderRequest,
    resources: session.resources,
    webFiles: session.webFiles || {},
    legacyFiles: session.files,
    storedConfig: session.config,
    fileBindings: session.cliInvocation?.fileBindings,
    repairInvalidComparisonHeight: Number(session.version) <= 33
  });
  const migratedConfig = migratePersistedGalleryConfig(
    mergedGuiConfig(session, projection)
  );
  const projectedFiles = {
    ...projection.files,
    linearSeqs: (projection.files.linearSeqs || []).map((record) => ({ ...record })),
    linearComparisons: (projection.files.linearComparisons || []).map((comparison) => ({
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
  const state = buildStateFacade(session, projection, config);
  restoreConservationFiles(session, filesData, state.circularConservation);
  const comparisonPlanSnapshot = projection.mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: state.linearComparisonPlan,
        sequences: filesData.linearSeqs,
        layout: state.linearRecordLayoutEnabled.value ? state.linearRecordRows : [],
        losatProgram: state.losatProgram.value,
        blastpMode: state.losat?.blastp?.mode
      })
    : null;
  const promotedCore = buildCanonicalRenderRequest({
    state,
    filesData,
    comparisonPlanSnapshot
  });
  const promoted = {
    ...session,
    config: cloneJson(config),
    renderRequest: promotedCore.renderRequest,
    resources: promotedCore.resources,
    webFiles: {
      ...(session.webFiles || {}),
      ...promotedCore.webFiles
    }
  };
  hydrateRecordPresentations(promoted.renderRequest, args);
  preserveComparisonResources(session, promoted);
  return promoted;
};

export const promoteGallerySessionToCurrent = (session) => {
  if (!isPlainObject(session) || !isPlainObject(session.renderRequest)) {
    throw new Error('Gallery session must contain a canonical renderRequest.');
  }
  if (!isPlainObject(session.resources)) {
    throw new Error('Gallery session must contain canonical resources.');
  }
  const schema = Number(session.renderRequest.schema);
  if (schema >= 5 && schema <= CANONICAL_REQUEST_SCHEMA) {
    if (Number(session.version) < 40) {
      const args = sessionArgs(session);
      const cliAuthored = isPlainObject(session?.config?.cliOptions)
        && !isPlainObject(session?.config?.linearRecordLayout);
      return cliAuthored
        ? promoteCliAuthoredSession(session, args)
        : promoteGuiAuthoredSession(session, args);
    }
    const promoted = cloneJson(session);
    if (isPlainObject(promoted.config)) {
      promoted.config = migratePersistedGalleryConfig(promoted.config);
    }
    return promoted;
  }
  if (![1, 2].includes(schema)) {
    throw new Error(`Unsupported canonical renderRequest schema: ${schema}.`);
  }
  const args = sessionArgs(session);
  const cliAuthored = isPlainObject(session?.config?.cliOptions)
    && !isPlainObject(session?.config?.linearRecordLayout);
  return cliAuthored
    ? promoteCliAuthoredSession(session, args)
    : promoteGuiAuthoredSession(session, args);
};
