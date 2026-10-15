// @ts-check
import { interpolateColor, resolveColorToHex } from '../utils/color-utils.js';

/** @import { FileLike } from './session-resource-backing.js' */

/**
 * One comparison ring: the source file it reads and how it is drawn.
 * @typedef {{
 *   sourceKey: string, fileName: string, sourceIndex: number, label: string, color: string, losat_gencode: number
 * }} ConservationSeriesEntry
 */
/**
 * A series entry as a Session or an older draft holds it: every field is unvalidated.
 * @typedef {{
 *   sourceKey?: unknown, fileName?: unknown, sourceIndex?: unknown, label?: unknown, name?: unknown,
 *   color?: unknown, losat_gencode?: unknown
 * }} RawConservationSeriesEntry
 */
/** @typedef {{ series?: RawConservationSeriesEntry[] | null, labels?: unknown }} CircularConservationLike */
/** @typedef {{ file: FileLike, sourceKey: string, sourceIndex: number, fileName: string, defaultLabel: string }} ConservationSourceDescriptor */
/** @typedef {ConservationSourceDescriptor & Omit<ConservationSeriesEntry, 'sourceKey' | 'fileName' | 'sourceIndex'> & { orderIndex: number }} OrderedConservationSource */

export const CONSERVATION_SLOT_MANAGER = 'circular_conservation';

const CONSERVATION_SERIES_COLOR_BLEND_TO_WHITE = 0.18;
const CONSERVATION_SERIES_BASE_COLORS = [
  '#4e79a7',
  '#f28e2b',
  '#59a14f',
  '#e15759',
  '#76b7b2',
  '#edc948',
  '#b07aa1',
  '#ff9da7',
  '#9c755f',
  '#bab0ac'
];
const CONSERVATION_SERIES_COLORS = CONSERVATION_SERIES_BASE_COLORS.map((color) =>
  interpolateColor(color, '#ffffff', CONSERVATION_SERIES_COLOR_BLEND_TO_WHITE)
);

/**
 * @template T
 * @param {T | T[] | null | undefined} files
 * @returns {T[]}
 */
export const normalizeFileList = (files) => (Array.isArray(files) ? files.filter(Boolean) : (files ? [files] : []));

/**
 * @template {FileLike} T
 * @param {(T | null)[] | null | undefined} files An empty place is kept as one.
 * @param {CircularConservationLike | null | undefined} circularConservation
 * @returns {(T | null)[]}
 */
export const orderedOptionalConservationFiles = (files, circularConservation) => {
  const sourceFiles = Array.isArray(files) ? files : [];
  const series = Array.isArray(circularConservation?.series)
    ? circularConservation.series
    : [];
  if (series.length === 0) return [...sourceFiles];

  /** @type {Set<number>} */
  const used = new Set();
  /** @type {(T | null)[]} */
  const ordered = [];
  series.forEach((entry, seriesIndex) => {
    const fileName = String(entry?.fileName || '').trim();
    let index = fileName
      ? sourceFiles.findIndex(
          (file, candidateIndex) => (
            !used.has(candidateIndex) &&
            String(file?.name || '').trim() === fileName
          )
        )
      : -1;
    if (index < 0) {
      const sourceIndex = Number(entry?.sourceIndex);
      index = Number.isInteger(sourceIndex) ? sourceIndex : seriesIndex;
    }
    if (index >= 0 && index < sourceFiles.length && !used.has(index)) {
      used.add(index);
      ordered.push(sourceFiles[index] || null);
      return;
    }
    ordered.push(null);
  });
  sourceFiles.forEach((file, index) => {
    if (!used.has(index)) ordered.push(file || null);
  });
  return ordered;
};

/**
 * The one removal of a comparison series with its files: Remove in the
 * Circular comparison panel, and Load for a series whose BLAST table it could
 * not read (OV-893, E-OD30). Each list holds one place per series, in series
 * order; the place `index` goes from every list, and the series after it keep
 * their own files and take the place before (`sourceIndex`).
 * @template F
 * @template {object} S
 * @param {{ blasts: F[], fastas: F[], sequenceSources: F[], series: S[] }} lists
 * @param {number} index
 */
export const withoutConservationSeries = ({ blasts, fastas, sequenceSources, series }, index) => {
  /** @template T @param {T[]} list */
  const without = (list) => list.filter((_, place) => place !== index);
  const kept = without(series).map((entry, sourceIndex) => ({ ...entry, sourceIndex }));
  return {
    blasts: without(blasts),
    fastas: without(fastas),
    sequenceSources: without(sequenceSources),
    series: kept,
    labels: kept.map((entry) => /** @type {{ label?: unknown }} */ (entry).label).join(',')
  };
};

/**
 * @param {{ name?: unknown } | null | undefined} file A File, or a series entry's file name; only `name` is read.
 * @param {unknown} index
 */
export const defaultConservationSeriesLabel = (file, index) => {
  const raw = String(file?.name || `Series ${Number(index) + 1}`).trim();
  const withoutExtension = raw.replace(/\.[^.]+$/, '').trim();
  return withoutExtension || `Series ${Number(index) + 1}`;
};

/**
 * @param {unknown} value
 * @returns {string[]}
 */
export const parseConservationLabelText = (value) =>
  String(value || '')
    .split(/[\n,]+/)
    .map((label) => label.trim())
    .filter(Boolean);

/**
 * @param {unknown} value
 * @param {unknown} index
 * @returns {string}
 */
export const normalizeConservationSeriesColor = (value, index) => {
  const fallback = CONSERVATION_SERIES_COLORS[Number(index) % CONSERVATION_SERIES_COLORS.length];
  const resolved = resolveColorToHex(String(value || fallback).trim());
  const color = String(resolved || fallback).trim();
  const shortMatch = color.match(/^#([0-9a-fA-F]{3})$/);
  if (shortMatch) {
    return `#${shortMatch[1].split('').map((char) => char + char).join('').toLowerCase()}`;
  }
  return /^#[0-9a-fA-F]{6}$/.test(color) ? color.toLowerCase() : fallback;
};

/**
 * @param {unknown} value
 * @param {number} [fallback]
 */
const normalizePositiveInteger = (value, fallback = 1) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : fallback;
};

/**
 * @param {FileLike} file `normalizeFileList` drops falsy entries
 * @param {unknown} index
 */
const fileSignature = (file, index) => {
  const fileName = String(file?.name || `source_${Number(index) + 1}`).trim();
  const size = Number.isFinite(Number(file?.size)) ? Number(file.size) : 0;
  const lastModified = Number.isFinite(Number(file?.lastModified)) ? Number(file.lastModified) : 0;
  return `${fileName}|${size}|${lastModified}`;
};

/**
 * @param {FileLike | FileLike[] | null | undefined} sourceFiles
 * @returns {ConservationSourceDescriptor[]}
 */
export const conservationSourceDescriptors = (sourceFiles) => {
  /** @type {Map<string, number>} */
  const seen = new Map();
  return normalizeFileList(sourceFiles).map((file, index) => {
    const signature = fileSignature(file, index);
    const duplicateIndex = seen.get(signature) || 0;
    seen.set(signature, duplicateIndex + 1);
    const fileName = String(file?.name || `source_${Number(index) + 1}`).trim();
    const sourceKey = `${signature}|${duplicateIndex}`;
    return {
      file,
      sourceKey,
      sourceIndex: index,
      fileName,
      defaultLabel: defaultConservationSeriesLabel(file, index)
    };
  });
};

/**
 * @param {ConservationSourceDescriptor} descriptor
 * @param {RawConservationSeriesEntry | null} previousEntry
 * @param {number} orderIndex
 * @param {string | null} [legacyLabel]
 * @returns {ConservationSeriesEntry}
 */
const mergeSeriesEntry = (descriptor, previousEntry, orderIndex, legacyLabel = null) => {
  const rawLabel = previousEntry?.label ?? previousEntry?.name ?? legacyLabel ?? descriptor.defaultLabel;
  return {
    sourceKey: descriptor.sourceKey,
    fileName: descriptor.fileName,
    sourceIndex: descriptor.sourceIndex,
    label: String(rawLabel || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
    color: normalizeConservationSeriesColor(previousEntry?.color, orderIndex),
    losat_gencode: normalizePositiveInteger(previousEntry?.losat_gencode, 1)
  };
};

/**
 * @param {{
 *   sourceFiles: FileLike | FileLike[] | null | undefined,
 *   previousSeries: RawConservationSeriesEntry[] | null | undefined,
 *   legacyLabels: string[] | null | undefined
 * }} options
 * @returns {ConservationSeriesEntry[]}
 */
export const reconcileConservationSeries = ({
  sourceFiles,
  previousSeries,
  legacyLabels
}) => {
  const descriptors = conservationSourceDescriptors(sourceFiles);
  const descriptorByKey = new Map(descriptors.map((descriptor) => [descriptor.sourceKey, descriptor]));
  /** @type {Map<string, ConservationSourceDescriptor[]>} */
  const descriptorsByFileName = new Map();
  descriptors.forEach((descriptor) => {
    if (!descriptorsByFileName.has(descriptor.fileName)) descriptorsByFileName.set(descriptor.fileName, []);
    // Set just above when absent.
    /** @type {ConservationSourceDescriptor[]} */ (descriptorsByFileName.get(descriptor.fileName)).push(descriptor);
  });

  /** @type {Set<string>} */
  const usedKeys = new Set();
  /** @type {ConservationSeriesEntry[]} */
  const nextSeries = [];
  const previous = Array.isArray(previousSeries) ? previousSeries : [];
  const fallbackLabels = Array.isArray(legacyLabels) ? legacyLabels : [];

  previous.forEach((entry) => {
    if (!entry || typeof entry !== 'object') return;
    const sourceKey = String(entry.sourceKey || '').trim();
    let descriptor = sourceKey ? descriptorByKey.get(sourceKey) : null;
    if (!descriptor) {
      const fileName = String(entry.fileName || '').trim();
      const candidates = descriptorsByFileName.get(fileName) || [];
      descriptor = candidates.find((candidate) => !usedKeys.has(candidate.sourceKey)) || null;
    }
    if (!descriptor || usedKeys.has(descriptor.sourceKey)) return;
    usedKeys.add(descriptor.sourceKey);
    nextSeries.push(mergeSeriesEntry(
      descriptor,
      entry,
      nextSeries.length,
      fallbackLabels[descriptor.sourceIndex]
    ));
  });

  descriptors.forEach((descriptor) => {
    if (usedKeys.has(descriptor.sourceKey)) return;
    usedKeys.add(descriptor.sourceKey);
    nextSeries.push(mergeSeriesEntry(
      descriptor,
      null,
      nextSeries.length,
      fallbackLabels[descriptor.sourceIndex]
    ));
  });

  return nextSeries;
};

// D12: a row added from a GenBank or DDBJ file takes the label the file names
// itself (first record's DEFINITION, then organism), read by the Python ring
// reader. Only the row's file-name default is replaced, so a typed label wins;
// a FASTA file names no label and keeps the file-name default.
/**
 * @param {{
 *   series: ConservationSeriesEntry[] | null | undefined,
 *   sourceFiles: FileLike | FileLike[] | null | undefined,
 *   file: FileLike | null | undefined,
 *   recordLabel: unknown
 * }} options
 */
export const applyComparisonSequenceRecordLabel = ({ series, sourceFiles, file, recordLabel }) => {
  const label = String(recordLabel ?? '').trim();
  if (!label || !file || !Array.isArray(series)) return false;
  const descriptor = conservationSourceDescriptors(sourceFiles).find((candidate) => candidate.file === file);
  const entry = descriptor
    ? series.find((candidate) => candidate?.sourceKey === descriptor.sourceKey)
    : null;
  if (!entry || !descriptor || String(entry.label ?? '').trim() !== descriptor.defaultLabel) return false;
  entry.label = label;
  return true;
};

/**
 * @param {FileLike | FileLike[] | null | undefined} sourceFiles
 * @param {CircularConservationLike | null | undefined} circularConservation
 * @returns {OrderedConservationSource[]}
 */
export const orderedConservationSources = (sourceFiles, circularConservation) => {
  const descriptors = conservationSourceDescriptors(sourceFiles);
  const legacyLabels = parseConservationLabelText(circularConservation?.labels);
  const descriptorByKey = new Map(descriptors.map((descriptor) => [descriptor.sourceKey, descriptor]));
  /** @type {Map<string, ConservationSourceDescriptor[]>} */
  const descriptorsByFileName = new Map();
  descriptors.forEach((descriptor) => {
    if (!descriptorsByFileName.has(descriptor.fileName)) descriptorsByFileName.set(descriptor.fileName, []);
    // Set just above when absent.
    /** @type {ConservationSourceDescriptor[]} */ (descriptorsByFileName.get(descriptor.fileName)).push(descriptor);
  });

  /** @type {Set<string>} */
  const usedKeys = new Set();
  /** @type {OrderedConservationSource[]} */
  const entries = [];
  const series = Array.isArray(circularConservation?.series) ? circularConservation.series : [];
  series.forEach((entry, seriesIndex) => {
    if (!entry || typeof entry !== 'object') return;
    const sourceKey = String(entry.sourceKey || '').trim();
    let descriptor = sourceKey ? descriptorByKey.get(sourceKey) : null;
    if (!descriptor) {
      const fileName = String(entry.fileName || '').trim();
      const candidates = descriptorsByFileName.get(fileName) || [];
      descriptor = candidates.find((candidate) => !usedKeys.has(candidate.sourceKey)) || null;
    }
    if (!descriptor) {
      const sourceIndex = Number(entry.sourceIndex);
      const fallbackIndex = Number.isInteger(sourceIndex) ? sourceIndex : seriesIndex;
      const indexed = descriptors[fallbackIndex];
      descriptor = indexed && !usedKeys.has(indexed.sourceKey) ? indexed : null;
    }
    if (!descriptor || usedKeys.has(descriptor.sourceKey)) return;
    usedKeys.add(descriptor.sourceKey);
    entries.push({
      ...descriptor,
      label: String(entry.label || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
      color: normalizeConservationSeriesColor(entry.color, entries.length),
      losat_gencode: normalizePositiveInteger(entry.losat_gencode, 1),
      orderIndex: entries.length
    });
  });

  descriptors.forEach((descriptor) => {
    if (usedKeys.has(descriptor.sourceKey)) return;
    usedKeys.add(descriptor.sourceKey);
    const fallbackLabel = legacyLabels[descriptor.sourceIndex] || descriptor.defaultLabel;
    entries.push({
      ...descriptor,
      label: String(fallbackLabel || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
      color: normalizeConservationSeriesColor(null, entries.length),
      losat_gencode: 1,
      orderIndex: entries.length
    });
  });

  return entries;
};

/**
 * @template T
 * @param {T[]} series
 * @param {unknown} index
 * @param {unknown} direction
 */
export const moveConservationSeriesEntry = (series, index, direction) => {
  const idx = Number(index);
  const step = Math.sign(Number(direction));
  if (!Array.isArray(series) || !Number.isInteger(idx) || step === 0) return false;
  const target = idx + step;
  if (idx < 0 || idx >= series.length || target < 0 || target >= series.length) return false;
  const [entry] = series.splice(idx, 1);
  series.splice(target, 0, entry);
  return true;
};

/**
 * @param {{ label?: unknown, fileName?: unknown } | null | undefined} entry
 * @param {unknown} index
 * @param {Set<string>} [existingIds]
 */
export const safeConservationSlotId = (entry, index, existingIds = new Set()) => {
  const source = String(entry?.label || entry?.fileName || `comparison_${Number(index) + 1}`).trim();
  const cleaned = source
    .replace(/\.[^.]+$/, '')
    .replace(/[^A-Za-z0-9_]+/g, '_')
    .replace(/^_+|_+$/g, '')
    .slice(0, 40) || `comparison_${Number(index) + 1}`;
  let id = `conservation_${cleaned}`;
  let suffix = 2;
  while (existingIds.has(id)) {
    id = `conservation_${cleaned}_${suffix}`;
    suffix += 1;
  }
  existingIds.add(id);
  return id;
};

/** @param {{ renderer?: unknown, params?: { managed?: unknown } | null } | null | undefined} slot */
export const isManagedConservationSlot = (slot) => (
  slot?.renderer === 'sequence_conservation' &&
  String(slot?.params?.managed || '').trim() === CONSERVATION_SLOT_MANAGER
);
