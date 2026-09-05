import {
  buildCircularTrackSlotSpec,
  parseCircularTrackSlotSpec
} from './circular-track-slots.js';
import {
  buildLinearTrackSlotSpec,
  parseLinearTrackSlotSpec
} from './linear-track-slots.js';
import { encodeAnnotationTable } from './annotations/table-codec.js';
import { base64ToBytes } from '../services/byte-utils.js';

const SAFE_SHELL_TOKEN_RE = /^[A-Za-z0-9_@%+=:,./-]+$/;

const SOURCE_PROVENANCE_UNAVAILABLE =
  'Source recipe unavailable: this session does not contain sufficient source-file provenance.';

class SourceRecipeUnavailable extends Error {}

const normalizePath = (value) => String(value ?? '').trim();

const isPlainObject = (value) => Boolean(
  value && typeof value === 'object' && !Array.isArray(value)
);

const visibleBasename = (value) => {
  const name = String(value || '')
    .replace(/\\/g, '/')
    .split('/')
    .pop()
    .replace(/[\u0000-\u001f\u007f]/g, '')
    .trim();
  return name && name !== '.' && name !== '..' ? name : '';
};

const fallbackNameFromPath = (path) => {
  const name = normalizePath(path).split('/').filter(Boolean).pop();
  return name || 'file';
};

const normalizeMetadataEntries = (fileMetadata) => {
  if (!fileMetadata) return [];
  if (fileMetadata instanceof Map) return Array.from(fileMetadata.entries());
  if (Array.isArray(fileMetadata)) {
    return fileMetadata
      .map((entry) => {
        if (Array.isArray(entry) && entry.length >= 2) return [entry[0], entry[1]];
        if (entry && typeof entry === 'object') return [entry.path, entry];
        return null;
      })
      .filter(Boolean);
  }
  if (typeof fileMetadata === 'object') return Object.entries(fileMetadata);
  return [];
};

const normalizeFileMetadata = (fileMetadata) => {
  const map = new Map();
  normalizeMetadataEntries(fileMetadata).forEach(([pathRaw, entryRaw]) => {
    const path = normalizePath(pathRaw || entryRaw?.path);
    if (!path) return;
    const entry = entryRaw && typeof entryRaw === 'object' ? entryRaw : {};
    const name = String(entry.name || fallbackNameFromPath(path)).trim() || fallbackNameFromPath(path);
    const kind = String(entry.kind || '').trim() === 'generated' ? 'generated' : 'uploaded';
    const slot = String(entry.slot || '').trim();
    map.set(path, {
      path,
      name,
      kind,
      slot
    });
  });
  return map;
};

const filenameCollisionKey = (name) => String(name || '')
  .replace(/[A-Z]/g, (letter) => letter.toLowerCase());

const suffixedFilename = (name, suffix) => {
  const dot = name.lastIndexOf('.');
  return dot > 0
    ? `${name.slice(0, dot)}-${suffix}${name.slice(dot)}`
    : `${name}-${suffix}`;
};

const allocateVisibleFilenames = ({ sourceMetadata, exactMetadata, sourceAvailable }) => {
  const source = new Map(Array.from(sourceMetadata, ([path, entry]) => [path, { ...entry }]));
  const exact = new Map(Array.from(exactMetadata, ([path, entry]) => [path, { ...entry }]));
  const used = new Map();
  const reserveFixed = (namespace, path, entry) => {
    const name = visibleBasename(entry.name);
    if (!name) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: an uploaded input has no safe visible filename.'
      );
    }
    const key = filenameCollisionKey(name);
    const owner = `${namespace}:${path}`;
    if (used.has(key) && used.get(key) !== owner) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: multiple uploaded inputs use the visible filename "${name}".`
      );
    }
    entry.name = name;
    used.set(key, owner);
  };
  const reserveGenerated = (namespace, path, entry, index) => {
    const base = visibleBasename(entry.name) || visibleBasename(path) || `helper-${index}`;
    let name = base;
    let suffix = 2;
    while (used.has(filenameCollisionKey(name))) {
      name = suffixedFilename(base, suffix);
      suffix += 1;
    }
    entry.name = name;
    used.set(filenameCollisionKey(name), `${namespace}:${path}`);
  };

  const namespaces = sourceAvailable ? [['source', source], ['exact', exact]] : [['exact', exact]];
  namespaces.forEach(([namespace, metadata]) => metadata.forEach((entry, path) => {
    if (entry.kind !== 'generated') reserveFixed(namespace, path, entry);
  }));
  let generatedIndex = 0;
  namespaces.forEach(([namespace, metadata]) => metadata.forEach((entry, path) => {
    if (entry.kind !== 'generated') return;
    generatedIndex += 1;
    reserveGenerated(namespace, path, entry, generatedIndex);
  }));
  return { source, exact };
};

const formatHelperFileList = (files) => {
  const names = Array.from(new Set((Array.isArray(files) ? files : [])
    .map((file) => String(file?.name || fallbackNameFromPath(file?.path)).trim())
    .filter(Boolean)));
  if (names.length === 0) return 'generated helper files';
  if (names.length <= 3) return names.join(', ');
  return `${names.slice(0, 3).join(', ')}, and ${names.length - 3} more`;
};

const selectorText = (selector) => {
  if (!isPlainObject(selector)) return '';
  if (selector.kind === 'recordId') return String(selector.value || '').trim();
  if (selector.kind === 'recordIndex') {
    const index = Number(selector.index);
    return Number.isInteger(index) && index >= 0 ? `#${index + 1}` : '';
  }
  return '';
};

const tsvCell = (value) => {
  const text = String(value ?? '');
  return /[\t\r\n"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
};

const tsv = (columns, rows) => `${[
  columns.join('\t'),
  ...rows.map((row) => columns.map((column) => tsvCell(row[column])).join('\t'))
].join('\n')}\n`;

const collectBindingNameHints = (value, hints) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => collectBindingNameHints(entry, hints));
    return;
  }
  if (!isPlainObject(value)) return;
  const resourceId = String(value.resourceId || '').trim();
  const name = visibleBasename(value.name);
  if (resourceId && name && !hints.has(resourceId)) hints.set(resourceId, name);
  Object.values(value).forEach((entry) => collectBindingNameHints(entry, hints));
};

const originalResourceNameHints = (webFiles) => {
  const hints = new Map();
  const stored = isPlainObject(webFiles?.resourceOriginalNames)
    ? webFiles.resourceOriginalNames
    : {};
  Object.entries(stored).forEach(([resourceId, name]) => {
    const visible = visibleBasename(name);
    if (resourceId && visible) hints.set(resourceId, visible);
  });
  collectBindingNameHints(webFiles?.bindings, hints);
  return hints;
};

const normalizeNameHints = (hints) => new Map(
  hints instanceof Map
    ? hints
    : (isPlainObject(hints) ? Object.entries(hints) : [])
);

const resourceBytes = (descriptor) => {
  if (!isPlainObject(descriptor)) return null;
  if (descriptor.encoding === 'base64' && typeof descriptor.data === 'string') {
    return base64ToBytes(descriptor.data);
  }
  if (descriptor.encoding === 'utf8' && typeof descriptor.data === 'string') {
    return new TextEncoder().encode(descriptor.data);
  }
  return null;
};

const createRecipeFiles = (resources, webFiles, generatedFileNameHints) => {
  const descriptors = isPlainObject(resources) ? resources : {};
  const originalNames = originalResourceNameHints(webFiles);
  const generatedNames = normalizeNameHints(generatedFileNameHints);
  const fileMetadata = new Map();
  const generatedFiles = [];
  const pathsByResourceId = new Map();
  const visibleNameOwners = new Map();
  const recordCountCache = new Map();
  let generatedIndex = 0;

  const reserveName = (name, owner) => {
    const collisionKey = name.replace(/[A-Z]/g, (letter) => letter.toLowerCase());
    const existing = visibleNameOwners.get(collisionKey);
    if (existing && existing !== owner) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: multiple required files use the visible filename "${name}".`
      );
    }
    visibleNameOwners.set(collisionKey, owner);
  };

  const resourcePath = (resourceIdRaw, {
    source = false,
    fallback = 'helper.tsv',
    slot = '',
    nameHintSlot = ''
  } = {}) => {
    const resourceId = String(resourceIdRaw || '').trim();
    const descriptor = descriptors[resourceId];
    if (!resourceId || !isPlainObject(descriptor)) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: canonical resource "${resourceId || '(missing)'}" is unavailable.`
      );
    }
    if (pathsByResourceId.has(resourceId)) return pathsByResourceId.get(resourceId);

    const originalName = originalNames.get(resourceId) || '';
    if (source && !originalName) throw new SourceRecipeUnavailable(SOURCE_PROVENANCE_UNAVAILABLE);
    const generated = !originalName;
    const name = originalName || visibleBasename(
      generatedNames.get(resourceId) || generatedNames.get(nameHintSlot) || fallback
    );
    if (!name) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: a generated helper has no trustworthy visible filename.'
      );
    }
    if (!generated) reserveName(name, resourceId);
    const path = `/run-info-resource-${resourceId}`;
    const resolvedSlot = generated
      ? (slot || `generatedFiles.source_recipe.${resourceId}`)
      : `resources.${resourceId}`;
    fileMetadata.set(path, {
      path,
      name,
      slot: resolvedSlot,
      kind: generated ? 'generated' : 'uploaded',
      resourceId
    });
    if (generated) {
      const data = resourceBytes(descriptor);
      if (!data) {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: generated helper file "${name}" cannot be materialized.`
        );
      }
      generatedFiles.push({ path, name, slot: resolvedSlot, data });
    }
    pathsByResourceId.set(resourceId, path);
    return path;
  };

  const generatedTextPath = (nameRaw, data, slot) => {
    generatedIndex += 1;
    const name = visibleBasename(nameRaw) || `helper-${generatedIndex}.tsv`;
    const path = `/run-info-generated-${generatedIndex}`;
    const resolvedSlot = slot || `generatedFiles.source_recipe.helper_${generatedIndex}`;
    fileMetadata.set(path, { path, name, slot: resolvedSlot, kind: 'generated' });
    generatedFiles.push({
      path,
      name,
      slot: resolvedSlot,
      ...(typeof data === 'function'
        ? { buildData: data }
        : { data: String(data ?? '') })
    });
    return path;
  };

  const visibleName = (path) => fileMetadata.get(path)?.name || fallbackNameFromPath(path);
  const resourceRecordCount = (resourceIdRaw, kind) => {
    const resourceId = String(resourceIdRaw || '').trim();
    const cacheKey = `${resourceId}:${kind}`;
    if (recordCountCache.has(cacheKey)) return recordCountCache.get(cacheKey);
    const bytes = resourceBytes(descriptors[resourceId]);
    let count = 0;
    if (bytes) {
      const text = new TextDecoder().decode(bytes);
      count = kind === 'genbank'
        ? (text.match(/^LOCUS\s+/gm) || []).length
        : (text.match(/^>/gm) || []).length;
    }
    recordCountCache.set(cacheKey, count);
    return count;
  };
  return {
    resourcePath,
    generatedTextPath,
    visibleName,
    resourceRecordCount,
    fileMetadata,
    generatedFiles
  };
};

const DIAGRAM_OPTION_FIELDS = new Set([
  'configOverrides', 'tracks', 'output', 'colors', 'selectedFeaturesSet', 'featureShapes',
  'dinucleotide', 'window', 'step', 'depthWindow', 'depthStep', 'depthTracks',
  'plotTitle', 'plotTitleFontSize', 'evalue', 'bitscore', 'identity', 'alignmentLength',
  'featureVisibilityTableFile', 'labelWhitelistFile', 'qualifierPriorityFile',
  'qualifierPriorityTable', 'labelOverrideFile', 'annotations', 'pairwiseMatchStyle',
  'keepFullDefinitionWithPlotTitle', 'species', 'strain', 'conservationBlastFiles',
  'conservationFastaFiles', 'conservationReference', 'conservationLabels',
  'conservationColors', 'conservationRingWidth', 'conservationRingGap'
]);

const coverObject = (coverage, value, path, projected, metadata = []) => {
  if (!isPlainObject(value)) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: semantic field "${path}" is malformed.`
    );
  }
  const projectedKeys = new Set(projected);
  const metadataKeys = new Set(metadata);
  Object.keys(value).forEach((key) => {
    const fieldPath = path ? `${path}.${key}` : key;
    if (projectedKeys.has(key)) coverage.consumedPaths.add(fieldPath);
    else if (metadataKeys.has(key)) coverage.metadataPaths.add(fieldPath);
    else {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: semantic field "${fieldPath}" has no current CLI projection.`
      );
    }
  });
};

const coverResourceRef = (coverage, value, path) => {
  if (value === null || value === undefined) return;
  coverObject(coverage, value, path, ['resourceId', 'representation']);
};

const coverTrackScalar = (coverage, value, path) => {
  if (isPlainObject(value)) coverObject(coverage, value, path, ['value', 'unit']);
};

const coverTrackSlots = (coverage, tracks, mode) => {
  coverObject(coverage, tracks, 'diagramOptions.tracks', [
    'circularTrackSlots', 'circularTrackAxisIndex', 'linearTrackSlots',
    'linearTrackAxisIndex', 'centerReservedRadius'
  ]);
  const slots = mode === 'circular' ? tracks.circularTrackSlots : tracks.linearTrackSlots;
  if (!Array.isArray(slots)) return;
  slots.forEach((slot, index) => {
    if (typeof slot === 'string') return;
    const path = `diagramOptions.tracks.${mode}TrackSlots[${index}]`;
    const scalarFields = mode === 'circular' ? ['radius', 'width'] : ['height', 'spacing'];
    coverObject(coverage, slot, path, [
      'kind', 'id', 'renderer', 'enabled', 'side', ...scalarFields,
      ...(mode === 'circular' ? ['innerGapPx', 'outerGapPx'] : []),
      'z', 'params'
    ]);
    scalarFields.forEach((field) => coverTrackScalar(coverage, slot[field], `${path}.${field}`));
    const params = slot.params;
    if (!isPlainObject(params)) return;
    coverObject(coverage, params, `${path}.params`, [
      'track_index', 'nt', 'positive_color', 'negative_color', 'tick_label_layout',
      'preset', 'lane_direction', 'source_index', 'set_id', 'marks', 'lane_gap_px',
      'padding_px', 'overflow', 'show_labels', 'anchor_slot', 'layer', 'cover_anchor',
      'legend_label'
    ]);
  });
};

const validateSchema6SemanticCoverage = (request) => {
  const coverage = { schema: 6, consumedPaths: new Set(), metadataPaths: new Set() };
  coverObject(coverage, request, '', [
    'mode', 'grouping', 'records', 'diagramOptions', 'layout', 'comparisons', 'output'
  ], ['schema']);
  const supportedGroupings = request.mode === 'circular' ? ['single', 'grid'] : ['single'];
  if (!supportedGroupings.includes(request.grouping) || Array.isArray(request.output)) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: batch or multiple-output semantics have no single current CLI projection.'
    );
  }

  (Array.isArray(request.records) ? request.records : []).forEach((record, index) => {
    const path = `records[${index}]`;
    coverObject(coverage, record, path, [
      'cardinality', 'source', 'selector', 'region', 'presentation'
    ], ['recordKey']);
    if (!['all', 'exactly_one'].includes(record.cardinality)) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: record cardinality "${String(record.cardinality || 'missing')}" has no current CLI projection.`
      );
    }
    coverObject(
      coverage,
      record.source,
      `${path}.source`,
      record.source?.kind === 'gffFasta'
        ? ['kind', 'gffResourceId', 'fastaResourceId']
        : ['kind', 'resourceId']
    );
    if (record.selector) {
      if (!['recordId', 'recordIndex'].includes(record.selector.kind)) {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: record selector "${String(record.selector.kind || 'unknown')}" has no current CLI projection.`
        );
      }
      coverObject(
        coverage,
        record.selector,
        `${path}.selector`,
        record.selector.kind === 'recordIndex' ? ['kind', 'index'] : ['kind', 'value']
      );
    }
    if (record.region) {
      coverObject(coverage, record.region, `${path}.region`, [
        'selector', 'start', 'end', 'reverseComplement'
      ]);
      if (record.region.selector) {
        if (!['recordId', 'recordIndex'].includes(record.region.selector.kind)) {
          throw new SourceRecipeUnavailable(
            `Source recipe unavailable: region selector "${String(record.region.selector.kind || 'unknown')}" has no current CLI projection.`
          );
        }
        coverObject(
          coverage,
          record.region.selector,
          `${path}.region.selector`,
          record.region.selector.kind === 'recordIndex' ? ['kind', 'index'] : ['kind', 'value']
        );
      }
    }
    coverObject(coverage, record.presentation, `${path}.presentation`, [
      'label', 'subtitle', 'reverseComplement', 'gridRow', 'gridColumn'
    ]);
    ['gridRow', 'gridColumn'].forEach((field) => {
      const value = record.presentation[field];
      if (value !== null && value !== undefined && value !== '' && gridCoordinate(value) === '') {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: record placement "${path}.presentation.${field}" is invalid.`
        );
      }
    });
  });

  const options = request.diagramOptions;
  coverObject(coverage, options, 'diagramOptions', DIAGRAM_OPTION_FIELDS);
  coverObject(coverage, options.configOverrides || {}, 'diagramOptions.configOverrides', Object.keys(
    options.configOverrides || {}
  ));
  coverTrackSlots(coverage, options.tracks || {}, request.mode);
  coverObject(coverage, options.output || {}, 'diagramOptions.output', ['legend', 'plotTitlePosition']);
  coverObject(coverage, options.colors || {}, 'diagramOptions.colors', [
    'colorTable', 'colorTableFile', 'defaultColors', 'defaultColorsPalette', 'defaultColorsFile'
  ]);
  ['colorTable', 'colorTableFile', 'defaultColors', 'defaultColorsFile'].forEach((field) => (
    coverResourceRef(coverage, options.colors?.[field], `diagramOptions.colors.${field}`)
  ));
  coverObject(coverage, options.featureShapes || {}, 'diagramOptions.featureShapes', Object.keys(
    options.featureShapes || {}
  ));
  [
    'featureVisibilityTableFile', 'labelWhitelistFile', 'qualifierPriorityFile',
    'qualifierPriorityTable', 'labelOverrideFile'
  ].forEach((field) => coverResourceRef(coverage, options[field], `diagramOptions.${field}`));
  (Array.isArray(options.depthTracks) ? options.depthTracks : []).forEach((track, index) => {
    const path = `diagramOptions.depthTracks[${index}]`;
    coverObject(coverage, track, path, [
      'source', 'label', 'color', 'height', 'largeTickInterval', 'smallTickInterval',
      'tickFontSize'
    ]);
    const refs = Array.isArray(track.source) ? track.source : [track.source];
    refs.forEach((ref, refIndex) => coverResourceRef(
      coverage, ref, `${path}.source${Array.isArray(track.source) ? `[${refIndex}]` : ''}`
    ));
  });
  ['conservationBlastFiles', 'conservationFastaFiles'].forEach((field) => {
    (Array.isArray(options[field]) ? options[field] : []).forEach((ref, index) => (
      coverResourceRef(coverage, ref, `diagramOptions.${field}[${index}]`)
    ));
  });
  if (options.annotations !== undefined) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: canonical annotation semantics have no fully lossless current CLI projection.'
    );
  }

  const layoutFields = request.mode === 'linear'
    ? ['recordGapPx']
    : [
        'multiRecordSizeMode', 'multiRecordMinRadiusRatio', 'multiRecordColumnGapRatio',
        'multiRecordRowGapRatio', 'multiRecordPositions'
      ];
  coverObject(coverage, request.layout || {}, 'layout', layoutFields);
  (Array.isArray(request.comparisons) ? request.comparisons : []).forEach((comparison, index) => {
    const path = `comparisons[${index}]`;
    if (comparison?.kind !== 'nucleotideBlast') {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: comparison "${String(comparison?.kind || 'unknown')}" has no lossless current CLI projection.`
      );
    }
    coverObject(coverage, comparison, path, [
      'kind', 'resourceId', 'queryRecordIndex', 'subjectRecordIndex'
    ]);
  });
  coverObject(coverage, request.output, 'output', [
    'prefix', 'formats', 'overwrite', 'interactiveMetadataPolicy'
  ]);
  if (request.output?.interactiveMetadataPolicy !== 'auto') {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: the interactive metadata policy has no current CLI projection.'
    );
  }
  return {
    schema: coverage.schema,
    consumedPaths: Array.from(coverage.consumedPaths).sort(),
    metadataPaths: Array.from(coverage.metadataPaths).sort()
  };
};

const sourceSpec = (record) => {
  const source = isPlainObject(record?.source) ? record.source : {};
  if (source.kind === 'genbank') {
    return { kind: 'genbank', ids: [String(source.resourceId || '').trim()] };
  }
  if (source.kind === 'gffFasta') {
    return {
      kind: 'gffFasta',
      ids: [
        String(source.gffResourceId || '').trim(),
        String(source.fastaResourceId || '').trim()
      ]
    };
  }
  throw new SourceRecipeUnavailable(
    `Source recipe unavailable: source kind "${String(source.kind || 'unknown')}" has no CLI recipe.`
  );
};

const circularLayoutRows = (request, records) => {
  const positions = Array.isArray(request.layout?.multiRecordPositions)
    ? request.layout.multiRecordPositions
    : [];
  if (positions.length === 0) return [];
  const rows = Array(records.length).fill(null);
  positions.forEach((position) => {
    const text = String(position || '').trim();
    const separator = text.lastIndexOf('@');
    const selector = separator > 0 ? text.slice(0, separator) : '';
    const row = Number(separator > 0 ? text.slice(separator + 1) : NaN);
    let recordIndex = -1;
    const indexMatch = selector.match(/^#(\d+)$/);
    if (indexMatch) {
      const index = Number(indexMatch[1]) - 1;
      if (Number.isInteger(index) && index >= 0 && index < records.length) recordIndex = index;
    } else {
      const matches = records
        .map((record, index) => (
          [selectorText(record?.selector), selectorText(record?.region?.selector)].includes(selector)
            ? index
            : -1
        ))
        .filter((index) => index >= 0);
      if (matches.length === 1) [recordIndex] = matches;
    }
    if (
      recordIndex < 0 ||
      !Number.isInteger(row) ||
      row < 1 ||
      rows[recordIndex] !== null
    ) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: circular record placement cannot be projected into a records table.'
      );
    }
    rows[recordIndex] = row;
  });
  if (rows.some((row) => row === null)) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: circular record placement is incomplete for a records table.'
    );
  }
  return rows;
};

const gridCoordinate = (value) => {
  if (value === null || value === undefined || value === '') return '';
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : '';
};

const appendInputArgs = (args, request, files) => {
  const records = Array.isArray(request.records) ? request.records : [];
  if (records.length === 0) {
    throw new SourceRecipeUnavailable('Source recipe unavailable: the committed render has no source records.');
  }
  const specs = records.map(sourceSpec);
  if (new Set(specs.map((spec) => spec.kind)).size !== 1) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: the committed render mixes source kinds that the CLI cannot combine.'
    );
  }
  const inputKind = specs[0].kind;
  const sourcePaths = specs.map((spec, recordIndex) => spec.ids.map((resourceId, partIndex) => (
    files.resourcePath(resourceId, {
      source: true,
      fallback: inputKind === 'genbank'
        ? `record-${recordIndex + 1}.gbk`
        : (partIndex === 0 ? `record-${recordIndex + 1}.gff3` : `record-${recordIndex + 1}.fasta`)
    })
  )));
  const sourceKeys = sourcePaths.map((paths) => paths.join('\0'));
  const hasDuplicateSources = new Set(sourceKeys).size !== sourceKeys.length;
  const recordNeedsTable = records.some((record) => (
    gridCoordinate(record.presentation?.gridColumn) !== ''
    || (request.mode === 'circular' && (
      Boolean(record.region)
      || Boolean(record.presentation?.label)
      || Boolean(record.presentation?.subtitle)
      || Boolean(record.presentation?.reverseComplement)
      || gridCoordinate(record.presentation?.gridRow) !== ''
    ))
  ));
  const circularSelectorNeedsTable = request.mode === 'circular' && records.some((record, index) => (
    Boolean(record.region?.selector || record.selector)
    && files.resourceRecordCount(
      specs[index].ids.at(-1),
      inputKind === 'genbank' ? 'genbank' : 'fasta'
    ) !== 1
  ));
  const recordsTableRequired = hasDuplicateSources || circularSelectorNeedsTable || recordNeedsTable;
  if (
    request.schema === 6
    && recordsTableRequired
    && records.some((record) => record.cardinality !== 'exactly_one')
  ) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: this record cardinality cannot be preserved by a CLI records table.'
    );
  }

  if (request.schema === 6 && !recordsTableRequired) {
    const loadComparison = request.mode === 'linear'
      && Array.isArray(request.comparisons)
      && request.comparisons.length > 0;
    const effectiveCardinality = request.mode === 'circular'
      ? 'all'
      : ((inputKind === 'gffFasta' || records.length > 1) && loadComparison ? 'first' : 'all');
    records.forEach((record, index) => {
      const canonicalCardinality = String(record.cardinality || '');
      const cliCardinality = (record.region?.selector || record.selector)
        ? 'exactly_one'
        : effectiveCardinality;
      if (
        canonicalCardinality !== cliCardinality
        && files.resourceRecordCount(
          specs[index].ids.at(-1),
          inputKind === 'genbank' ? 'genbank' : 'fasta'
        ) !== 1
      ) {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: record cardinality "${canonicalCardinality || 'missing'}" cannot be preserved by the current CLI.`
        );
      }
    });
  }

  if (recordsTableRequired) {
    const layoutRows = request.mode === 'circular'
      ? circularLayoutRows(request, records)
      : [];
    const columns = inputKind === 'genbank'
      ? ['gbk', 'record_label', 'record_subtitle', 'record_id', 'region', 'reverse_complement', 'order', 'row', 'column']
      : ['gff', 'fasta', 'record_label', 'record_subtitle', 'record_id', 'region', 'reverse_complement', 'order', 'row', 'column'];
    const rows = records.map((record, index) => {
      const region = isPlainObject(record.region) ? record.region : null;
      const selector = selectorText(region?.selector || record.selector);
      const row = {
        gbk: inputKind === 'genbank' ? files.visibleName(sourcePaths[index][0]) : '',
        gff: inputKind === 'gffFasta' ? files.visibleName(sourcePaths[index][0]) : '',
        fasta: inputKind === 'gffFasta' ? files.visibleName(sourcePaths[index][1]) : '',
        record_label: String(record.presentation?.label || ''),
        record_subtitle: String(record.presentation?.subtitle || ''),
        record_id: selector,
        region: region
          ? `${Number(region.start)}-${Number(region.end)}${region.reverseComplement ? ':rc' : ''}`
          : '',
        reverse_complement: region ? '0' : (record.presentation?.reverseComplement ? '1' : '0'),
        order: index + 1,
        row: gridCoordinate(record.presentation?.gridRow) || layoutRows[index] || '',
        column: gridCoordinate(record.presentation?.gridColumn)
      };
      return row;
    });
    const tablePath = files.generatedTextPath(
      'records.tsv',
      tsv(columns, rows),
      'generatedFiles.source_recipe.records_table'
    );
    args.push('--records_table', tablePath);
    return true;
  }

  if (inputKind === 'genbank') {
    args.push('--gbk', ...sourcePaths.map(([gbk]) => gbk));
  } else {
    args.push('--gff', ...sourcePaths.map(([gff]) => gff));
    args.push('--fasta', ...sourcePaths.map(([, fasta]) => fasta));
  }

  if (request.mode !== 'linear') return false;
  const selectors = records.map((record) => selectorText(record.selector));
  if (selectors.some(Boolean)) args.push(...selectors.flatMap((value) => ['--record_id', value]));
  const labels = records.map((record) => String(record.presentation?.label || ''));
  if (labels.some(Boolean)) args.push(...labels.flatMap((value) => ['--record_label', value]));
  const subtitles = records.map((record) => String(record.presentation?.subtitle || ''));
  if (subtitles.some(Boolean)) args.push(...subtitles.flatMap((value) => ['--record_subtitle', value]));
  const reverseFlags = records.map((record) => Boolean(record.presentation?.reverseComplement));
  if (reverseFlags.some(Boolean)) {
    args.push(...reverseFlags.flatMap((value) => ['--reverse_complement', value ? '1' : '0']));
  }
  records.forEach((record, index) => {
    if (!isPlainObject(record.region)) return;
    const selector = selectorText(record.region.selector) || `#${index + 1}`;
    args.push(
      '--region',
      `${selector}:${Number(record.region.start)}-${Number(record.region.end)}` +
        `${record.region.reverseComplement ? ':rc' : ''}`
    );
  });
  const rows = records.map((record) => Number(record.presentation?.gridRow));
  if (rows.some((row) => Number.isInteger(row) && row > 0)) {
    rows.forEach((row, index) => {
      if (Number.isInteger(row) && row > 0) args.push('--multi_record_position', `#${index + 1}@${row}`);
    });
  }
  return false;
};

const appendOption = (args, option, value) => {
  if (value === null || value === undefined || value === '') return;
  args.push(option, String(value));
};

const appendBooleanOption = (args, value, enabledOption, disabledOption = '') => {
  if (value === true && enabledOption) args.push(enabledOption);
  else if (value === false && disabledOption) args.push(disabledOption);
};

const CONFIG_VALUE_OPTIONS = Object.freeze({
  'objects.features.arrow_geometry.head_length_ratio': '--arrow_head_length_ratio',
  'objects.features.arrow_geometry.shaft_width_ratio': '--arrow_shaft_width_ratio',
  'objects.features.block_stroke_color': '--block_stroke_color',
  'objects.features.line_stroke_color': '--line_stroke_color',
  'labels.rendering': '--label_rendering',
  'objects.gc_content.mode': '--gc_content_mode',
  'objects.gc_content.min_percent': '--gc_content_min_percent',
  'objects.gc_content.max_percent': '--gc_content_max_percent',
  'objects.gc_content.large_tick_interval': '--gc_content_large_tick_interval',
  'objects.gc_content.small_tick_interval': '--gc_content_small_tick_interval',
  'objects.gc_content.tick_font_size': '--gc_content_tick_font_size',
  'objects.depth.fill_color': '--depth_color',
  'objects.depth.min_depth': '--depth_min',
  'objects.depth.max_depth': '--depth_max',
  'objects.depth.large_tick_interval': '--depth_large_tick_interval',
  'objects.depth.small_tick_interval': '--depth_small_tick_interval',
  'objects.depth.tick_font_size': '--depth_tick_font_size',
  'objects.scale.interval': '--scale_interval'
});

const CIRCULAR_CONFIG_VALUE_OPTIONS = Object.freeze({
  'objects.axis.circular.stroke_color': '--axis_stroke_color',
  'objects.definition.circular.font_size': '--definition_font_size',
  'canvas.circular.track_type': '--track_type',
  'labels.spacing.circular': '--circular_label_spacing',
  'labels.circular.placement': '--label_placement',
  'objects.ticks.tick_labels.font_size': '--tick_label_font_size',
  'labels.unified_adjustment.outer_labels.x_radius_offset': '--outer_label_x_radius_offset',
  'labels.unified_adjustment.outer_labels.y_radius_offset': '--outer_label_y_radius_offset',
  'labels.unified_adjustment.inner_labels.x_radius_offset': '--inner_label_x_radius_offset',
  'labels.unified_adjustment.inner_labels.y_radius_offset': '--inner_label_y_radius_offset'
});

const LINEAR_CONFIG_VALUE_OPTIONS = Object.freeze({
  'objects.axis.linear.stroke_color': '--axis_stroke_color',
  'labels.spacing.linear': '--linear_label_spacing',
  'labels.linear.placement': '--label_placement',
  'labels.linear.rotation': '--label_rotation',
  'canvas.linear.track_layout': '--track_layout',
  'canvas.linear.track_axis_gap': '--track_axis_gap',
  'canvas.linear.comparison_height': '--comparison_height',
  'canvas.linear.default_gc_height': '--gc_height',
  'canvas.linear.depth_height': '--depth_height',
  'objects.scale.style': '--scale_style',
  'objects.scale.stroke_color': '--scale_stroke_color',
  'objects.scale.label_color': '--ruler_label_color',
  'objects.scale.stroke_width': '--scale_stroke_width',
  'objects.blast_match.style': '--pairwise_match_style'
});

const SHARED_LENGTH_OPTIONS = Object.freeze({
  'objects.features.block_stroke_width': '--block_stroke_width',
  'objects.features.line_stroke_width': '--line_stroke_width',
  'objects.legends.color_rect_size': '--legend_box_size',
  'objects.legends.font_size': '--legend_font_size'
});

const CIRCULAR_LENGTH_OPTIONS = Object.freeze({
  'objects.axis.circular.stroke_width': '--axis_stroke_width',
  'labels.font_size': '--label_font_size'
});

const LINEAR_LENGTH_OPTIONS = Object.freeze({
  'objects.axis.linear.stroke_width': '--axis_stroke_width',
  'objects.definition.linear.font_size': '--definition_font_size',
  'canvas.linear.default_cds_height': '--feature_height',
  'objects.scale.font_size': '--scale_font_size',
  'objects.scale.ruler_label_font_size': '--ruler_label_font_size',
  'labels.font_size.linear': '--label_font_size'
});

const appendConfigOverrides = (args, request) => {
  const options = request.diagramOptions || {};
  const overrides = isPlainObject(options.configOverrides) ? options.configOverrides : {};
  const consumed = new Set();
  const take = (path) => {
    if (!Object.hasOwn(overrides, path)) return undefined;
    consumed.add(path);
    return overrides[path];
  };

  Object.entries({
    ...CONFIG_VALUE_OPTIONS,
    ...(request.mode === 'linear' ? LINEAR_CONFIG_VALUE_OPTIONS : CIRCULAR_CONFIG_VALUE_OPTIONS)
  }).forEach(([path, option]) => appendOption(args, option, take(path)));

  appendBooleanOption(args, take('canvas.show_gc'), '--gc', '--no-gc');
  appendBooleanOption(args, take('canvas.show_skew'), '--skew', '--no-skew');
  take('canvas.show_depth');
  appendBooleanOption(args, take('canvas.strandedness'), '--separate_strands');
  appendBooleanOption(args, take('canvas.resolve_overlaps'), '--resolve_overlaps');
  appendBooleanOption(
    args, take('objects.gc_content.show_axis'),
    '--show_gc_content_axis', '--hide_gc_content_axis'
  );
  appendBooleanOption(
    args, take('objects.gc_content.show_ticks'),
    '--show_gc_content_ticks', '--hide_gc_content_ticks'
  );
  appendBooleanOption(args, take('objects.depth.normalize'), '--depth_log_scale', '--no_depth_log_scale');
  appendBooleanOption(args, take('objects.depth.show_axis'), '--show_depth_axis', '--hide_depth_axis');
  appendBooleanOption(args, take('objects.depth.show_ticks'), '--show_depth_ticks', '--hide_depth_ticks');
  appendBooleanOption(args, take('objects.depth.share_axis'), '--share_depth_axis');
  appendBooleanOption(args, take('objects.scale.show'), '', '--hide_scale');

  const blacklist = take('labels.filtering.blacklist_keywords');
  if (Array.isArray(blacklist) && blacklist.length > 0) {
    appendOption(args, '--label_blacklist', blacklist.join(','));
  } else if (blacklist !== undefined && !Array.isArray(blacklist)) {
    throw new SourceRecipeUnavailable('Source recipe unavailable: the label blacklist is malformed.');
  }

  const scopePath = request.mode === 'linear' ? 'labels.linear.scope' : 'labels.circular.scope';
  const scope = take(scopePath);
  if (scope !== undefined) {
    if (request.mode === 'linear') {
      appendOption(args, '--show_labels', scope);
    } else {
      const value = { outer: 'out', both: 'both', none: 'none' }[String(scope)];
      if (!value) throw new SourceRecipeUnavailable(`Source recipe unavailable: unsupported label scope "${scope}".`);
      appendOption(args, '--labels', value);
    }
  }

  if (request.mode === 'linear') {
    appendBooleanOption(args, take('objects.definition.linear.show_replicon'), '--show_replicon');
    appendBooleanOption(args, take('objects.definition.linear.show_accession'), '', '--hide_accession');
    appendBooleanOption(args, take('objects.definition.linear.show_length'), '', '--hide_length');
    appendBooleanOption(args, take('canvas.linear.align_center'), '--align_center');
    appendBooleanOption(args, take('canvas.linear.keep_definition_left_aligned'), '--keep_definition_left_aligned');
    appendBooleanOption(args, take('canvas.linear.ruler_on_axis'), '--ruler_on_axis');
    appendBooleanOption(args, take('canvas.linear.normalize_length'), '--normalize_length');
  }

  const lengthOptions = {
    ...SHARED_LENGTH_OPTIONS,
    ...(request.mode === 'linear' ? LINEAR_LENGTH_OPTIONS : CIRCULAR_LENGTH_OPTIONS)
  };
  Object.entries(lengthOptions).forEach(([prefix, option]) => {
    const shortValue = take(`${prefix}.short`);
    const longValue = take(`${prefix}.long`);
    if (shortValue === undefined && longValue === undefined) return;
    if (shortValue === undefined || longValue === undefined || Number(shortValue) !== Number(longValue)) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: size-specific setting "${prefix}" cannot be expressed by the current CLI.`
      );
    }
    appendOption(args, option, shortValue);
  });

  if (request.mode === 'circular') {
    const definitionInterval = take('objects.definition.circular.interval');
    if (definitionInterval !== undefined) {
      const fontSize = overrides['objects.definition.circular.font_size'];
      if (fontSize === undefined || Number(definitionInterval) !== Math.trunc(Number(fontSize) + 2)) {
        throw new SourceRecipeUnavailable(
          'Source recipe unavailable: the custom circular definition interval has no CLI option.'
        );
      }
    }
    const configuredTitleFontSize = take('objects.definition.circular.plot_title_font_size');
    if (
      configuredTitleFontSize !== undefined
      && Number(configuredTitleFontSize) !== Number(options.plotTitleFontSize)
    ) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: conflicting plot-title font settings cannot be expressed by the current CLI.'
      );
    }
  } else {
    ['name', 'subtitle', 'replicon', 'accession', 'length'].forEach((kind) => {
      const prefix = `objects.definition.linear.line_styles.${kind}`;
      const fields = [
        ['font_size', 'size'],
        ['font_weight', 'weight'],
        ['fill', 'color']
      ];
      const values = fields
        .map(([field, key]) => [key, take(`${prefix}.${field}`)])
        .filter(([, value]) => value !== undefined);
      if (values.length > 0) {
        appendOption(
          args,
          '--definition_line_style',
          `${kind}:${values.map(([key, value]) => `${key}=${value}`).join(',')}`
        );
      }
    });
  }

  const unsupported = Object.keys(overrides).filter((path) => !consumed.has(path));
  if (unsupported.length > 0) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: setting "${unsupported[0]}" has no current CLI projection.`
    );
  }
};

const referencedResourceId = (value) => String(value?.resourceId || '').trim();

const appendResourceOptions = (args, options, files) => {
  const colors = isPlainObject(options.colors) ? options.colors : {};
  appendOption(args, '-p', colors.defaultColorsPalette);
  const resourceOptions = [
    ['-t', colors.colorTable || colors.colorTableFile, 'specific-colors.tsv'],
    ['-d', colors.defaultColors || colors.defaultColorsFile, 'default-colors.tsv'],
    ['--feature_visibility_table', options.featureVisibilityTableFile, 'feature-visibility.tsv'],
    ['--label_whitelist', options.labelWhitelistFile, 'label-whitelist.tsv'],
    ['--qualifier_priority', options.qualifierPriorityFile || options.qualifierPriorityTable, 'qualifier-priority.tsv'],
    ['--label_table', options.labelOverrideFile, 'label-overrides.tsv']
  ];
  resourceOptions.forEach(([option, ref, fallback]) => {
    const resourceId = referencedResourceId(ref);
    if (resourceId) appendOption(args, option, files.resourcePath(resourceId, { fallback }));
  });
  if (Array.isArray(options.annotations?.sets) && options.annotations.sets.length > 0) {
    appendOption(
      args,
      '--annotation_table',
      files.generatedTextPath(
        'annotations.tsv',
        encodeAnnotationTable(options.annotations.sets),
        'generatedFiles.source_recipe.annotation_table'
      )
    );
  }
};

const appendDepthOptions = (args, options, files) => {
  const tracks = Array.isArray(options.depthTracks) ? options.depthTracks : [];
  tracks.forEach((track, trackIndex) => {
    const refs = Array.isArray(track?.source) ? track.source : [track?.source];
    const paths = refs.map((ref, recordIndex) => {
      const resourceId = referencedResourceId(ref);
      return resourceId
        ? files.resourcePath(resourceId, {
            fallback: `depth-${trackIndex + 1}-record-${recordIndex + 1}.tsv`
          })
        : '';
    });
    if (!paths.some(Boolean)) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: depth track ${trackIndex + 1} has no source file.`
      );
    }
    args.push('--depth_track', ...paths);
  });
  const appendTrackValues = (option, field) => {
    if (tracks.length === 0) return;
    const values = tracks.map((track) => {
      const value = track?.[field];
      return value === null || value === undefined || value === '' ? 'auto' : String(value);
    });
    appendOption(args, option, values[0]);
    if (values.length > 1) args.push(...values.slice(1));
  };
  appendTrackValues('--depth_track_label', 'label');
  appendTrackValues('--depth_track_color', 'color');
  if (options.mode === 'linear') appendTrackValues('--depth_track_height', 'height');
  appendTrackValues('--depth_track_large_tick_interval', 'largeTickInterval');
  appendTrackValues('--depth_track_small_tick_interval', 'smallTickInterval');
  appendTrackValues('--depth_track_tick_font_size', 'tickFontSize');
};

const formatCanonicalTrackScalar = (value, { allowFactor = false } = {}) => {
  if (!isPlainObject(value)) return value;
  const numeric = Number(value.value);
  const unit = String(value.unit || '').trim().toLowerCase();
  if (!Number.isFinite(numeric) || (unit !== 'px' && !(allowFactor && unit === 'factor'))) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: a canonical track dimension cannot be represented by the current CLI.'
    );
  }
  return unit === 'px' ? `${numeric}px` : String(numeric);
};

const TRACK_PARAM_KEYS = Object.freeze({
  circular: {
    features: new Set(['lane_direction', 'legend_label']),
    ticks: new Set(['tick_label_layout', 'preset', 'legend_label']),
    dinucleotide_content: new Set(['nt', 'legend_label']),
    dinucleotide_skew: new Set(['nt', 'positive_color', 'negative_color', 'legend_label']),
    depth: new Set(['track_index', 'legend_label']),
    sequence_conservation: new Set(['track_index', 'source_index', 'legend_label']),
    annotations: new Set([
      'set_id', 'marks', 'lane_gap_px', 'padding_px', 'overflow', 'show_labels',
      'anchor_slot', 'layer', 'cover_anchor', 'legend_label'
    ]),
    spacer: new Set(['legend_label'])
  },
  linear: {
    features: new Set(['legend_label']),
    dinucleotide_content: new Set(['nt', 'legend_label']),
    dinucleotide_skew: new Set(['nt', 'positive_color', 'negative_color', 'legend_label']),
    depth: new Set(['track_index', 'legend_label']),
    annotations: new Set([
      'set_id', 'marks', 'lane_gap_px', 'padding_px', 'overflow', 'show_labels',
      'anchor_slot', 'layer', 'cover_anchor', 'legend_label'
    ]),
    spacer: new Set(['legend_label'])
  }
});

const validateTrackParams = (mode, renderer, params) => {
  if (!isPlainObject(params)) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: a canonical ${mode} track slot has malformed params.`
    );
  }
  const allowed = TRACK_PARAM_KEYS[mode]?.[renderer];
  if (!allowed) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: ${mode} track renderer "${renderer}" is not accepted by the current CLI.`
    );
  }
  const unsupported = Object.keys(params).find((key) => !allowed.has(key));
  if (unsupported) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: ${mode} track field "params.${unsupported}" cannot be projected losslessly.`
    );
  }
};

const validateStructuredTrackSlot = (slot, mode) => {
  if (!isPlainObject(slot)) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: a canonical ${mode} track slot is malformed.`
    );
  }
  if (slot.kind !== `${mode}TrackSlot`) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: a canonical ${mode} track slot has an invalid kind.`
    );
  }
  const id = typeof slot.id === 'string' ? slot.id.trim() : '';
  const renderer = typeof slot.renderer === 'string' ? slot.renderer.trim() : '';
  const allowedSides = mode === 'circular'
    ? new Set(['inside', 'outside', 'overlay'])
    : new Set(['above', 'below', 'overlay']);
  if (
    !id || !renderer || typeof slot.enabled !== 'boolean' || !Number.isInteger(slot.z)
    || (slot.side !== null && !allowedSides.has(slot.side))
  ) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: a canonical ${mode} track slot has invalid identity fields.`
    );
  }
  validateTrackParams(mode, renderer, slot.params);
  return slot;
};

const trackSlotForRecipe = (slot, mode, options) => {
  if (typeof slot === 'string') {
    const text = slot.trim();
    if (!text) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: a canonical ${mode} track slot cannot be empty.`
      );
    }
    try {
      const parsed = mode === 'circular'
        ? parseCircularTrackSlotSpec(text)
        : parseLinearTrackSlotSpec(text);
      const atIndex = text.indexOf('@');
      const canonicalHead = `${parsed.id}:${parsed.renderer}`;
      if ((atIndex < 0 ? text : text.slice(0, atIndex)).trim() !== canonicalHead) {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: a canonical ${mode} track slot would change identity in the current CLI projection.`
        );
      }
      if (atIndex >= 0 && text.slice(atIndex + 1).split(',').some((entry) => {
        const equalsIndex = entry.indexOf('=');
        return equalsIndex <= 0 || !entry.slice(equalsIndex + 1).trim();
      })) {
        throw new SourceRecipeUnavailable(
          `Source recipe unavailable: a canonical ${mode} track slot has a malformed option.`
        );
      }
      validateTrackParams(mode, parsed.renderer, parsed.params || {});
    } catch (error) {
      if (error instanceof SourceRecipeUnavailable) throw error;
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: a canonical ${mode} track slot is not accepted by the current CLI.`
      );
    }
    return text;
  }

  validateStructuredTrackSlot(slot, mode);
  try {
    if (mode === 'circular') {
      return buildCircularTrackSlotSpec(
        circularTrackSlotForRecipe(slot),
        options.dinucleotide,
        options.configOverrides?.['canvas.circular.track_type']
      );
    }
    const parsed = parseLinearTrackSlotSpec(slot);
    validateTrackParams('linear', parsed.renderer, parsed.params || {});
    return buildLinearTrackSlotSpec(linearTrackSlotForRecipe(parsed), { includeEnabled: true });
  } catch (error) {
    if (error instanceof SourceRecipeUnavailable) throw error;
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: a canonical ${mode} track slot cannot be projected losslessly.`
    );
  }
};

const circularTrackSlotForRecipe = (slot) => ({
  ...slot,
  radius: formatCanonicalTrackScalar(slot?.radius, { allowFactor: true }),
  width: formatCanonicalTrackScalar(slot?.width, { allowFactor: true })
});

const linearTrackSlotForRecipe = (slot) => ({
  ...slot,
  height: formatCanonicalTrackScalar(slot?.height),
  spacing: formatCanonicalTrackScalar(slot?.spacing)
});

const appendTrackOptions = (args, request) => {
  const options = request.diagramOptions || {};
  const tracks = isPlainObject(options.tracks) ? options.tracks : {};
  appendOption(args, '--center_reserved_radius', tracks.centerReservedRadius);
  const field = request.mode === 'circular' ? 'circularTrackSlots' : 'linearTrackSlots';
  const axisField = request.mode === 'circular' ? 'circularTrackAxisIndex' : 'linearTrackAxisIndex';
  const slots = tracks[field];
  if (slots === null || slots === undefined) {
    if (tracks[axisField] !== null && tracks[axisField] !== undefined) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: ${request.mode} track axis placement has no explicit track slots.`
      );
    }
    return;
  }
  if (!Array.isArray(slots)) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: canonical ${request.mode} track slots must be an array or null.`
    );
  }
  if (slots.length === 0) {
    throw new SourceRecipeUnavailable(
      `Source recipe unavailable: an explicit empty ${request.mode} track list means no tracks, which the current CLI cannot express.`
    );
  }
  const option = request.mode === 'circular' ? '--circular_track_slot' : '--linear_track_slot';
  slots.forEach((slot) => appendOption(args, option, trackSlotForRecipe(slot, request.mode, options)));
  appendOption(
    args,
    request.mode === 'circular' ? '--circular_track_axis_index' : '--linear_track_axis_index',
    tracks[axisField]
  );
};

const appendComparisonOptions = (args, request, files) => {
  const options = request.diagramOptions || {};
  if (request.mode === 'circular') {
    const blasts = Array.isArray(options.conservationBlastFiles)
      ? options.conservationBlastFiles
      : [];
    if (blasts.length === 0) return;
    const fastas = Array.isArray(options.conservationFastaFiles)
      ? options.conservationFastaFiles
      : [];
    const labels = Array.isArray(options.conservationLabels) ? options.conservationLabels : [];
    const colors = Array.isArray(options.conservationColors) ? options.conservationColors : [];
    args.push(
      '--conservation_blast',
      ...blasts.map((ref, index) => files.resourcePath(
        referencedResourceId(ref),
        {
          fallback: `conservation-${index + 1}.tsv`,
          slot: `generatedFiles.source_recipe.conservation_blasts[${index}]`,
          nameHintSlot: `generatedFiles.circular_conservation_blasts[${index}]`
        }
      ))
    );
    if (fastas.some(Boolean)) {
      args.push(
        '--conservation_fasta',
        ...fastas.map((ref, index) => {
        const resourceId = referencedResourceId(ref);
          return resourceId
            ? files.resourcePath(resourceId, { fallback: `comparison-${index + 1}.fasta` })
            : '';
        })
      );
    }
    if (labels.length > 0) args.push('--conservation_labels', ...labels.map(String));
    if (colors.length > 0) args.push('--conservation_colors', ...colors.map(String));
    appendOption(args, '--conservation_reference', options.conservationReference);
    appendOption(args, '--conservation_ring_width', options.conservationRingWidth);
    appendOption(args, '--conservation_ring_gap', options.conservationRingGap);
    return;
  }

  const comparisons = Array.isArray(request.comparisons) ? request.comparisons : [];
  if (comparisons.length === 0) return;
  const recordCount = Array.isArray(request.records) ? request.records.length : 0;
  if (comparisons.some((comparison) => (
    !Number.isInteger(comparison?.queryRecordIndex)
    || comparison.queryRecordIndex < 0
    || comparison.queryRecordIndex >= recordCount
    || !Number.isInteger(comparison?.subjectRecordIndex)
    || comparison.subjectRecordIndex < 0
    || comparison.subjectRecordIndex >= recordCount
  ))) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: a comparison cannot be bound losslessly to its source records.'
    );
  }
  if (comparisons.some((comparison) => comparison.kind !== 'nucleotideBlast')) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: this committed protein comparison has no lossless source CLI projection.'
    );
  }
  const drawable = comparisons;
  const adjacent = drawable.every((comparison, index) => (
    Number(comparison.queryRecordIndex) === index
    && Number(comparison.subjectRecordIndex) === index + 1
  ));
  if (adjacent) {
    const paths = drawable.map((comparison, index) => files.resourcePath(
      referencedResourceId(comparison),
      {
        fallback: `comparison-${index + 1}.tsv`,
        slot: `generatedFiles.source_recipe.comparison_blasts[${index}]`,
        nameHintSlot: `generatedFiles.losat_blasts[${index}]`
      }
    ));
    if (paths.length > 0) args.push('-b', ...paths);
    return;
  }
  const rows = drawable.map((comparison, index) => {
    const path = files.resourcePath(referencedResourceId(comparison), {
      fallback: `comparison-${index + 1}.tsv`,
      slot: `generatedFiles.source_recipe.comparison_blasts[${index}]`,
      nameHintSlot: `generatedFiles.losat_blasts[${index}]`
    });
    return {
      blastPath: path,
      query: `#${Number(comparison.queryRecordIndex) + 1}`,
      subject: `#${Number(comparison.subjectRecordIndex) + 1}`
    };
  });
  appendOption(
    args,
    '--comparisons_table',
    files.generatedTextPath(
      'comparisons.tsv',
      (finalNameForPath) => tsv(
        ['blast', 'query', 'subject'],
        rows.map(({ blastPath, ...row }) => ({
          ...row,
          blast: finalNameForPath(blastPath)
        }))
      ),
      'generatedFiles.source_recipe.comparisons_table'
    )
  );
};

const appendDiagramOptions = (args, request, files) => {
  const options = request.diagramOptions || {};
  if (Array.isArray(request.output) && request.output.length !== 1) {
    throw new SourceRecipeUnavailable(
      'Source recipe unavailable: per-record batch output names have no single current CLI projection.'
    );
  }
  const output = Array.isArray(request.output) ? request.output[0] : request.output;
  appendOption(args, '-o', output?.prefix);
  if (Array.isArray(output?.formats) && output.formats.length > 0) {
    appendOption(args, '-f', output.formats.join(','));
  }
  if (output?.overwrite === true) args.push('--overwrite');
  appendOption(args, '-l', options.output?.legend);
  if (request.mode === 'linear' || options.output?.plotTitlePosition !== 'none') {
    appendOption(args, '--plot_title_position', options.output?.plotTitlePosition);
  }
  if (options.selectedFeaturesSet !== null && options.selectedFeaturesSet !== undefined) {
    if (!Array.isArray(options.selectedFeaturesSet)) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: selected features must be an array or null.'
      );
    }
    if (!options.selectedFeaturesSet.every((feature) => (
      typeof feature === 'string' && feature.trim().length > 0
    ))) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: selected features must contain only non-empty strings.'
      );
    }
    if (options.selectedFeaturesSet.length > 0) {
      args.push('-k', options.selectedFeaturesSet.join(','));
    }
  }
  Object.entries(isPlainObject(options.featureShapes) ? options.featureShapes : {})
    .forEach(([feature, shape]) => appendOption(args, '--feature_shape', `${feature}=${shape}`));
  appendOption(args, '-n', options.dinucleotide);
  appendOption(args, '-w', options.window);
  appendOption(args, '-s', options.step);
  appendOption(args, '--depth_window', options.depthWindow);
  appendOption(args, '--depth_step', options.depthStep);
  appendOption(args, '--plot_title', options.plotTitle);
  appendOption(args, '--plot_title_font_size', options.plotTitleFontSize);
  ['evalue', 'bitscore', 'identity'].forEach((field) => appendOption(args, `--${field}`, options[field]));
  appendOption(args, '--alignment_length', options.alignmentLength);
  if (request.mode === 'circular') {
    appendOption(args, '--species', options.species);
    appendOption(args, '--strain', options.strain);
    appendBooleanOption(
      args,
      options.keepFullDefinitionWithPlotTitle,
      '--keep_full_definition_with_plot_title'
    );
  } else if (Object.hasOwn(options, 'pairwiseMatchStyle')) {
    const configured = options.configOverrides?.['objects.blast_match.style'];
    if (configured === undefined) appendOption(args, '--pairwise_match_style', options.pairwiseMatchStyle);
    else if (String(configured) !== String(options.pairwiseMatchStyle)) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: conflicting pairwise comparison styles cannot be expressed by the current CLI.'
      );
    }
  }
  appendResourceOptions(args, options, files);
  appendDepthOptions(args, { ...options, mode: request.mode }, files);
  appendTrackOptions(args, request);
  appendComparisonOptions(args, request, files);
};

const appendLayoutOptions = (args, request, recordsTableUsed = false) => {
  const layout = isPlainObject(request.layout) ? request.layout : {};
  if (request.mode === 'linear') {
    appendOption(args, '--linear_record_gap', layout.recordGapPx);
    return;
  }
  if (Object.keys(layout).length === 0) return;
  args.push('--multi_record_canvas');
  appendOption(args, '--multi_record_size_mode', layout.multiRecordSizeMode);
  appendOption(args, '--multi_record_min_radius_ratio', layout.multiRecordMinRadiusRatio);
  appendOption(args, '--multi_record_column_gap_ratio', layout.multiRecordColumnGapRatio);
  appendOption(args, '--multi_record_row_gap_ratio', layout.multiRecordRowGapRatio);
  if (!recordsTableUsed && Array.isArray(layout.multiRecordPositions)) {
    layout.multiRecordPositions.forEach((position) => appendOption(args, '--multi_record_position', position));
  }
};

export const buildSourceRecipe = ({
  renderRequest,
  resources,
  webFiles,
  generatedFileNameHints
} = {}) => {
  const mode = String(renderRequest?.mode || '').trim();
  const unavailable = (reason) => ({
    available: false,
    mode: mode === 'linear' ? 'linear' : 'circular',
    args: [],
    fileMetadata: new Map(),
    generatedFiles: [],
    semanticCoverage: null,
    unavailableReason: String(reason || 'Source recipe unavailable.')
  });
  if (!isPlainObject(renderRequest) || !['circular', 'linear'].includes(mode)) {
    return unavailable('Source recipe unavailable: the committed render request is missing or invalid.');
  }
  try {
    const semanticCoverage = renderRequest.schema === 6
      ? validateSchema6SemanticCoverage(renderRequest)
      : { schema: renderRequest.schema, consumedPaths: [], metadataPaths: [] };
    const files = createRecipeFiles(resources, webFiles, generatedFileNameHints);
    const args = [];
    const recordsTableUsed = appendInputArgs(args, renderRequest, files);
    appendDiagramOptions(args, renderRequest, files);
    appendConfigOverrides(args, renderRequest);
    appendLayoutOptions(args, renderRequest, recordsTableUsed);
    return {
      available: true,
      mode,
      args,
      fileMetadata: files.fileMetadata,
      generatedFiles: files.generatedFiles,
      semanticCoverage,
      unavailableReason: ''
    };
  } catch (error) {
    if (error instanceof SourceRecipeUnavailable) return unavailable(error.message);
    console.warn('Failed to project the committed render into a source CLI recipe.', error);
    return unavailable('Source recipe unavailable: the committed render cannot be represented safely by the current CLI.');
  }
};

export const quoteShellArg = (value) => {
  const token = String(value ?? '');
  if (token.length > 0 && SAFE_SHELL_TOKEN_RE.test(token)) return token;
  return `'${token.replace(/'/g, `'\\''`)}'`;
};

export const buildShellCommand = (tokens) => (Array.isArray(tokens) ? tokens : [])
  .map((token) => quoteShellArg(token))
  .join(' ');

export const formatElapsedMs = (elapsedMs) => {
  const ms = Number(elapsedMs);
  if (!Number.isFinite(ms) || ms < 0) return 'n/a';
  if (ms < 1000) return `${Math.round(ms)} ms`;
  if (ms < 60_000) return `${(ms / 1000).toFixed(ms < 10_000 ? 2 : 1)} s`;
  const minutes = Math.floor(ms / 60_000);
  const seconds = (ms - minutes * 60_000) / 1000;
  return `${minutes} min ${seconds.toFixed(1)} s`;
};

export const reproducibilityLabel = (level) => {
  const normalized = String(level || '').trim();
  if (normalized === 'source-unavailable') return 'Source recipe unavailable';
  if (normalized === 'canonical-request') return 'Exact typed request';
  if (normalized === 'exact-uploaded-files') return 'Source recipe ready';
  if (normalized === 'requires-helper-files') return 'Source recipe needs helper files';
  if (normalized === 'session-recommended') return 'Source recipe needs helper files';
  return 'Approximate command';
};

export const isCliInvocationSessionExportable = (invocation) => {
  if (!invocation || typeof invocation !== 'object') return false;
  if (invocation.sessionExportable === false) return false;
  const bindings = Array.isArray(invocation.fileBindings) ? invocation.fileBindings : [];
  return bindings.every((binding) => String(binding?.slot || '').startsWith('files.'));
};

const projectInvocation = ({ mode, args, fileMetadata, generatedBy }) => {
  const metadata = normalizeFileMetadata(fileMetadata);
  const helperFiles = [];
  const unresolvedFileArgs = [];
  const displayArgs = (Array.isArray(args) ? args : []).map((arg, argIndex) => {
    const token = String(arg ?? '');
    const meta = metadata.get(token);
    if (!meta) {
      if (token.startsWith('/')) unresolvedFileArgs.push({ argIndex, path: token });
      return token;
    }
    const displayName = meta.name || fallbackNameFromPath(token);
    if (meta.kind === 'generated') {
      helperFiles.push({ path: meta.path, name: displayName, slot: meta.slot || '' });
    }
    return displayName;
  });
  const bindingArgIndexesByName = new Map();
  displayArgs.forEach((token, index) => {
    if (!bindingArgIndexesByName.has(token)) bindingArgIndexesByName.set(token, []);
    bindingArgIndexesByName.get(token).push(index);
  });
  const fixedFileBindings = [];
  const seenBindingKeys = new Set();
  (Array.isArray(args) ? args : []).forEach((arg, originalIndex) => {
    const meta = metadata.get(String(arg ?? ''));
    if (!meta?.slot) return;
    const displayName = meta.name || fallbackNameFromPath(meta.path);
    const indexes = bindingArgIndexesByName.get(displayName) || [];
    const argIndex = indexes.shift();
    if (!Number.isInteger(argIndex)) return;
    const key = `${argIndex}:${meta.slot}`;
    if (seenBindingKeys.has(key)) return;
    seenBindingKeys.add(key);
    fixedFileBindings.push({ argIndex, slot: meta.slot, name: displayName, originalArgIndex: originalIndex });
  });
  const invocation = {
    schema: 1,
    mode,
    args: displayArgs,
    renderFormats: ['svg'],
    fileBindings: fixedFileBindings.map(({ originalArgIndex: _originalArgIndex, ...binding }) => binding),
    generatedBy,
    sessionExportable: fixedFileBindings.every((binding) => String(binding.slot).startsWith('files.'))
      && unresolvedFileArgs.length === 0
  };
  const commandArgs = ['gbdraw', mode, ...displayArgs];
  return {
    command: buildShellCommand(commandArgs),
    commandArgs,
    invocation,
    helperFiles,
    unresolvedFileArgs
  };
};

const finalizeGeneratedRecipeFiles = (generatedFiles, allocatedMetadata) => {
  const finalNameForPath = (pathRaw) => {
    const path = String(pathRaw || '').trim();
    const name = allocatedMetadata.get(path)?.name;
    if (!name) {
      throw new SourceRecipeUnavailable(
        'Source recipe unavailable: a generated helper dependency has no allocated filename.'
      );
    }
    return name;
  };
  return (Array.isArray(generatedFiles) ? generatedFiles : []).map((file) => {
    const path = String(file?.path || '').trim();
    const slot = String(file?.slot || '').trim();
    const name = finalNameForPath(path);
    const data = typeof file?.buildData === 'function'
      ? file.buildData(finalNameForPath)
      : file?.data;
    if (!(typeof data === 'string' || data instanceof Uint8Array)) {
      throw new SourceRecipeUnavailable(
        `Source recipe unavailable: generated helper file "${name}" cannot be materialized.`
      );
    }
    return { path, name, slot, data };
  });
};

export const buildRunInfo = ({
  mode,
  args,
  sourceRecipe,
  exactReplayArgs,
  fileMetadata,
  elapsedMs,
  resultCount,
  startedAtIso,
  generatedBy = 'gbdraw-web'
} = {}) => {
  const normalizedMode = String(mode || '').trim() === 'linear' ? 'linear' : 'circular';
  let sourceAvailable = sourceRecipe?.available !== false;
  let sourceUnavailableReason = String(sourceRecipe?.unavailableReason || 'Source recipe unavailable.');
  const sourceMetadata = normalizeFileMetadata(
    sourceRecipe?.fileMetadata || fileMetadata
  );
  const allExactMetadata = normalizeFileMetadata(fileMetadata);
  const exactPaths = new Set((Array.isArray(exactReplayArgs) ? exactReplayArgs : []).map(String));
  const exactMetadata = new Map(
    Array.from(allExactMetadata).filter(([path]) => exactPaths.has(path))
  );
  let allocatedMetadata;
  let finalizedGeneratedFiles = [];
  try {
    allocatedMetadata = allocateVisibleFilenames({
      sourceMetadata,
      exactMetadata,
      sourceAvailable
    });
    if (sourceAvailable) {
      finalizedGeneratedFiles = finalizeGeneratedRecipeFiles(
        sourceRecipe?.generatedFiles,
        allocatedMetadata.source
      );
      if (Array.isArray(sourceRecipe?.generatedFiles)) {
        sourceRecipe.generatedFiles.splice(
          0,
          sourceRecipe.generatedFiles.length,
          ...finalizedGeneratedFiles
        );
      }
    }
  } catch (error) {
    if (!(error instanceof SourceRecipeUnavailable)) throw error;
    sourceAvailable = false;
    sourceUnavailableReason = error.message;
    allocatedMetadata = allocateVisibleFilenames({
      sourceMetadata: new Map(),
      exactMetadata,
      sourceAvailable: false
    });
  }
  const source = sourceAvailable
    ? projectInvocation({
        mode: normalizedMode,
        args: Array.isArray(sourceRecipe?.args) ? sourceRecipe.args : args,
        fileMetadata: allocatedMetadata.source,
        generatedBy
      })
    : {
        command: '', commandArgs: [], helperFiles: [], unresolvedFileArgs: [], invocation: null
      };
  const exact = Array.isArray(exactReplayArgs)
    ? projectInvocation({
        mode: normalizedMode,
        args: exactReplayArgs,
        fileMetadata: allocatedMetadata.exact,
        generatedBy
      })
    : null;
  const recipeGeneratedHelpers = (Array.isArray(sourceRecipe?.generatedFiles)
    ? finalizedGeneratedFiles
    : [])
    .map(({ path, name, slot }) => ({ path, name, slot }))
    .filter((file) => file.path && file.name);
  const sourceHelpers = Array.from(new Map(
    [...source.helperFiles, ...recipeGeneratedHelpers]
      .map((helper) => [`${helper.path}:${helper.slot}`, helper])
  ).values());
  const helperFiles = Array.from(new Map(
    [...sourceHelpers, ...(exact?.helperFiles || [])]
      .map((helper) => [`${helper.path}:${helper.slot}`, helper])
  ).values());
  const hasLosatHelpers = sourceHelpers.some(
    (helper) => /losat|blast/i.test(`${helper.slot} ${helper.name}`)
  );
  const notes = [];
  let level = 'exact-uploaded-files';
  if (!sourceAvailable) {
    level = 'source-unavailable';
    notes.push(sourceUnavailableReason);
    if (exact) notes.push('Exact replay remains available from the committed canonical session.');
  } else if (sourceHelpers.length > 0) {
    level = hasLosatHelpers ? 'session-recommended' : 'requires-helper-files';
    notes.push(
      `The source recipe references browser-generated file(s): ${formatHelperFileList(sourceHelpers)}. ` +
      'Use "Download reproducibility files" to save those file(s) next to the uploaded inputs before running the command in a terminal.'
    );
    notes.push(
      'If you save a .gbdraw-session.json.gz and load it back in the web app, the uploaded inputs, results, and LOSAT cache are restored from the compressed JSON; these TSV files do not need to be saved separately for session restore.'
    );
  }
  if (source.unresolvedFileArgs.length > 0) {
    level = 'pseudo';
    notes.push('Some browser virtual paths could not be mapped to visible file names, so this command is only an approximation.');
  }
  if (notes.length === 0) {
    notes.push('The source recipe can be rerun with the uploaded file names shown here.');
  }
  if (exact?.helperFiles.length > 0) {
    notes.push(
      `Exact replay references ${formatHelperFileList(exact.helperFiles)}, which is included in "Download reproducibility files".`
    );
  }
  const bundleDescription = [
    sourceAvailable
      ? (
          sourceHelpers.length > 0
            ? `Source recipe helper files: ${formatHelperFileList(sourceHelpers)}.`
            : 'The Source recipe uses only the original source files.'
        )
      : 'The Source recipe is unavailable.',
    exact?.helperFiles.length > 0
      ? `Exact replay file: ${formatHelperFileList(exact.helperFiles)}.`
      : 'No Exact replay file is available.'
  ].join(' ');

  return {
    schema: 1,
    mode: normalizedMode,
    startedAtIso: startedAtIso || new Date().toISOString(),
    elapsedMs: Number.isFinite(Number(elapsedMs)) ? Number(elapsedMs) : 0,
    resultCount: Number.isFinite(Number(resultCount)) ? Number(resultCount) : 0,
    command: source.command,
    commandArgs: source.commandArgs,
    sessionCommand: exact?.command || '',
    invocation: source.invocation,
    sourceRecipe: {
      available: sourceAvailable,
      command: source.command,
      commandArgs: source.commandArgs,
      invocation: source.invocation,
      helperFiles: sourceHelpers,
      reproducibility: {
        status: !sourceAvailable
          ? 'unavailable'
          : (sourceHelpers.length > 0 ? 'requires-generated-helper-files' : 'ready-with-original-files')
      },
      unavailableReason: sourceAvailable
        ? ''
        : sourceUnavailableReason
    },
    exactReplay: {
      available: Boolean(exact),
      command: exact?.command || '',
      commandArgs: exact?.commandArgs || [],
      helperFiles: exact?.helperFiles || [],
      reproducibility: { status: exact ? 'available' : 'unavailable' }
    },
    helperFiles,
    downloadBundle: {
      available: helperFiles.length > 0,
      sourceRecipeFiles: sourceHelpers,
      exactReplayFiles: exact?.helperFiles || [],
      description: bundleDescription
    },
    reproducibility: {
      level,
      notes
    }
  };
};
