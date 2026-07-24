import { normalizePaletteColors } from '../app/color-utils.js';
import {
  buildCanonicalSessionRequest,
  projectCanonicalSessionRequest
} from './session-request.js';

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

const textToBase64 = (text) => {
  const bytes = new TextEncoder().encode(String(text));
  let binary = '';
  const chunkSize = 0x8000;
  for (let index = 0; index < bytes.length; index += chunkSize) {
    binary += String.fromCharCode(...bytes.subarray(index, index + chunkSize));
  }
  return btoa(binary);
};

const cachedLosatFile = (entry, fallbackName) => {
  const text = String(entry?.text || '');
  const bytes = new TextEncoder().encode(text);
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
};

const mergedGuiConfig = (session, projection) => {
  const saved = isPlainObject(session.config) ? session.config : {};
  const projected = isPlainObject(projection.config) ? projection.config : {};
  const ui = isPlainObject(session.ui) ? session.ui : {};
  const appliedPaletteColors = isPlainObject(ui.appliedPaletteColors)
    ? ui.appliedPaletteColors
    : null;
  return {
    ...projected,
    ...saved,
    form: { ...(projected.form || {}), ...(saved.form || {}) },
    adv: { ...(projected.adv || {}), ...(saved.adv || {}) },
    palette: String(ui.appliedPaletteName || saved.palette || projected.palette || 'default'),
    colors: appliedPaletteColors && Object.keys(appliedPaletteColors).length > 0
      ? cloneJson(appliedPaletteColors)
      : cloneJson(saved.colors || projected.colors || {}),
    annotationSets: Array.isArray(saved.annotationSets)
      ? cloneJson(saved.annotationSets)
      : cloneJson(projected.annotationSets || [])
  };
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
    blastSource: ref(config.blastSource || 'files'),
    losatProgram: ref(config.losatProgram || 'blastn'),
    losat: cloneJson(config.losat || { blastp: {} }),
    selectedOrthogroupAlignmentFeature: ref(
      session?.orthogroupState?.selectedOrthogroupAlignmentFeature || ''
    ),
    linearRecordLayoutEnabled: ref(Boolean(layout.enabled)),
    linearRecordGap: ref(layout.recordGap ?? 24),
    linearRecordRows: cloneJson(layout.rows || []),
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
  const promoted = {
    ...session,
    renderRequest: cloneJson(session.renderRequest),
    resources: session.resources
  };
  const renderRequest = promoted.renderRequest;
  renderRequest.schema = 3;
  renderRequest.diagramOptions = isPlainObject(renderRequest.diagramOptions)
    ? renderRequest.diagramOptions
    : {};
  const featureShapes = isPlainObject(renderRequest.diagramOptions.featureShapes)
    ? renderRequest.diagramOptions.featureShapes
    : {};
  const explicitRepeat = explicitRepeatRendering(args);
  featureShapes.repeat_region = explicitRepeat || 'underlay';
  renderRequest.diagramOptions.featureShapes = featureShapes;
  hydrateRecordPresentations(renderRequest, args);
  for (const comparison of renderRequest.comparisons || []) {
    const resourceId = String(comparison?.resourceId || '');
    if (resourceId && !promoted.resources[resourceId]) {
      throw new Error(`Canonical comparison resource is missing: ${resourceId}`);
    }
  }
  return promoted;
};

const promoteGuiAuthoredSession = (session, args) => {
  const projection = projectCanonicalSessionRequest({
    renderRequest: session.renderRequest,
    resources: session.resources,
    webFiles: session.webFiles || {},
    legacyFiles: session.files,
    fileBindings: session.cliInvocation?.fileBindings,
    repairInvalidComparisonHeight: Number(session.version) <= 33
  });
  const config = mergedGuiConfig(session, projection);
  const filesData = {
    ...projection.files,
    linearSeqs: (projection.files.linearSeqs || []).map((record) => ({ ...record })),
    linearComparisons: (projection.files.linearComparisons || []).map((comparison) => ({
      ...comparison
    }))
  };
  hydrateLinearFilePresentations(filesData, args);
  const state = buildStateFacade(session, projection, config);
  restoreConservationFiles(session, filesData, state.circularConservation);
  const promotedCore = buildCanonicalSessionRequest({ state, filesData });
  const promoted = {
    ...session,
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

export const promoteGallerySessionToCanonicalV3 = (session) => {
  if (!isPlainObject(session) || !isPlainObject(session.renderRequest)) {
    throw new Error('Gallery session must contain a canonical renderRequest.');
  }
  if (!isPlainObject(session.resources)) {
    throw new Error('Gallery session must contain canonical resources.');
  }
  const args = sessionArgs(session);
  const cliAuthored = isPlainObject(session?.config?.cliOptions);
  return cliAuthored
    ? promoteCliAuthoredSession(session, args)
    : promoteGuiAuthoredSession(session, args);
};
