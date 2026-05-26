import { resolveColorToHex } from './color-utils.js';

export const CONSERVATION_SLOT_MANAGER = 'circular_conservation';

const CONSERVATION_SERIES_COLORS = [
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

export const normalizeFileList = (files) => (Array.isArray(files) ? files.filter(Boolean) : (files ? [files] : []));

export const defaultConservationSeriesLabel = (file, index) => {
  const raw = String(file?.name || `Series ${Number(index) + 1}`).trim();
  const withoutExtension = raw.replace(/\.[^.]+$/, '').trim();
  return withoutExtension || `Series ${Number(index) + 1}`;
};

export const parseConservationLabelText = (value) =>
  String(value || '')
    .split(/[\n,]+/)
    .map((label) => label.trim())
    .filter(Boolean);

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

const fileSignature = (file, index) => {
  const fileName = String(file?.name || `source_${Number(index) + 1}`).trim();
  const size = Number.isFinite(Number(file?.size)) ? Number(file.size) : 0;
  const lastModified = Number.isFinite(Number(file?.lastModified)) ? Number(file.lastModified) : 0;
  return `${fileName}|${size}|${lastModified}`;
};

export const conservationSourceDescriptors = (sourceFiles) => {
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

const mergeSeriesEntry = (descriptor, previousEntry, orderIndex, legacyLabel = null) => {
  const rawLabel = previousEntry?.label ?? previousEntry?.name ?? legacyLabel ?? descriptor.defaultLabel;
  return {
    sourceKey: descriptor.sourceKey,
    fileName: descriptor.fileName,
    sourceIndex: descriptor.sourceIndex,
    label: String(rawLabel || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
    color: normalizeConservationSeriesColor(previousEntry?.color, orderIndex)
  };
};

export const reconcileConservationSeries = ({
  sourceFiles,
  previousSeries,
  legacyLabels
}) => {
  const descriptors = conservationSourceDescriptors(sourceFiles);
  const descriptorByKey = new Map(descriptors.map((descriptor) => [descriptor.sourceKey, descriptor]));
  const descriptorsByFileName = new Map();
  descriptors.forEach((descriptor) => {
    if (!descriptorsByFileName.has(descriptor.fileName)) descriptorsByFileName.set(descriptor.fileName, []);
    descriptorsByFileName.get(descriptor.fileName).push(descriptor);
  });

  const usedKeys = new Set();
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

export const orderedConservationSources = (sourceFiles, circularConservation) => {
  const descriptors = conservationSourceDescriptors(sourceFiles);
  const legacyLabels = parseConservationLabelText(circularConservation?.labels);
  const descriptorByKey = new Map(descriptors.map((descriptor) => [descriptor.sourceKey, descriptor]));
  const descriptorsByFileName = new Map();
  descriptors.forEach((descriptor) => {
    if (!descriptorsByFileName.has(descriptor.fileName)) descriptorsByFileName.set(descriptor.fileName, []);
    descriptorsByFileName.get(descriptor.fileName).push(descriptor);
  });

  const usedKeys = new Set();
  const entries = [];
  const series = Array.isArray(circularConservation?.series) ? circularConservation.series : [];
  series.forEach((entry) => {
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
    entries.push({
      ...descriptor,
      label: String(entry.label || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
      color: normalizeConservationSeriesColor(entry.color, entries.length),
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
      orderIndex: entries.length
    });
  });

  return entries;
};

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

export const isManagedConservationSlot = (slot) => (
  slot?.renderer === 'sequence_conservation' &&
  String(slot?.params?.managed || '').trim() === CONSERVATION_SLOT_MANAGER
);
