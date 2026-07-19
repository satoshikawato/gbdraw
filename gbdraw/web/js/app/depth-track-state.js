import { normalizeOptionalText } from './track-slot-display.js';

const DEPTH_TRACK_FALLBACK_COLORS = [
  '#4A90E2',
  '#E45756',
  '#2CA02C',
  '#F28E2B',
  '#9467BD',
  '#8C564B',
  '#17BECF',
  '#7F7F7F'
];

const normalizePositiveNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

const fallbackLabel = (index) => {
  const numericIndex = Math.max(0, Number(index) || 0);
  return numericIndex === 0 ? 'Depth' : `Depth ${numericIndex + 1}`;
};

const fallbackColor = (index, defaults = {}) => {
  if (typeof defaults.colorForIndex === 'function') {
    return String(defaults.colorForIndex(index) || '');
  }
  if (index === 0 && defaults.depthColor) return String(defaults.depthColor);
  return DEPTH_TRACK_FALLBACK_COLORS[index % DEPTH_TRACK_FALLBACK_COLORS.length];
};

export const depthFileBaseName = (file) => {
  const rawName = String(file?.name || '').split(/[\\/]/).pop().trim();
  if (!rawName) return '';
  return rawName.replace(/\.[^.]+$/, '').trim() || rawName;
};

export const depthLabelMatchesFile = (label, file) => {
  const labelText = String(label || '').toLowerCase();
  const fileText = depthFileBaseName(file).toLowerCase();
  if (!labelText || !fileText) return false;
  if (labelText.includes(fileText) || fileText.includes(labelText)) return true;
  const accessionTokens = fileText.match(/[a-z]{2,5}\d{5,}/g) || [];
  return accessionTokens.some((token) => labelText.includes(token));
};

export const depthFileSlotsFromValue = (value) => {
  if (Array.isArray(value)) return value.slice();
  return value ? [value] : [];
};

export const isRecordMajorDepthFileMatrix = (value) => (
  Array.isArray(value) && value.every((row) => Array.isArray(row))
);

export const parseDepthTrackIndexIdentity = (value, fieldName = 'track_index') => {
  let numeric = null;
  if (typeof value === 'number') {
    numeric = value;
  } else if (
    typeof value === 'string' &&
    /^[+-]?\d+$/.test(value.trim())
  ) {
    numeric = Number(value.trim());
  }
  if (!Number.isSafeInteger(numeric) || numeric < 0) {
    throw new Error(`${fieldName} must be a non-negative integer.`);
  }
  return numeric;
};

export const normalizeRecordMajorDepthFileRows = (value, recordCount = 1) => {
  const requestedCount = Math.max(1, Number(recordCount) || 1);
  const sourceIsMatrix = isRecordMajorDepthFileMatrix(value);
  const rows = sourceIsMatrix
    ? value.map((row) => depthFileSlotsFromValue(row))
    : Array.from(
        { length: requestedCount },
        () => depthFileSlotsFromValue(value)
      );
  const targetCount = Math.max(requestedCount, rows.length);
  while (rows.length < targetCount) rows.push([]);
  const width = depthTrackMatrixWidth(rows);
  return rows.map((row) => padDepthFileSlots(row, width));
};

export const representativeDepthFiles = (value) => {
  const rows = isRecordMajorDepthFileMatrix(value)
    ? normalizeRecordMajorDepthFileRows(value, value.length)
    : [depthFileSlotsFromValue(value)];
  const width = depthTrackMatrixWidth(rows);
  return Array.from({ length: width }, (_, trackIndex) => (
    rows.map((row) => row[trackIndex]).find(Boolean) || null
  ));
};

export const compactDepthFileSlots = (slots) => {
  const next = depthFileSlotsFromValue(slots);
  while (next.length > 0 && !next[next.length - 1]) {
    next.pop();
  }
  return next;
};

export const uploadedDepthFileCount = (value) => (
  depthFileSlotsFromValue(value).filter(Boolean).length
);

export const depthTrackMatrixWidth = (rows) => (
  (Array.isArray(rows) ? rows : []).reduce(
    (width, row) => Math.max(width, depthFileSlotsFromValue(row).length),
    0
  )
);

export const activeDepthTrackIndices = (rows) => {
  const active = new Set();
  (Array.isArray(rows) ? rows : []).forEach((row) => {
    depthFileSlotsFromValue(row).forEach((file, trackIndex) => {
      if (file) active.add(trackIndex);
    });
  });
  return Array.from(active).sort((left, right) => left - right);
};

export const depthTrackCoverageCount = (rows, trackIndex) => {
  const idx = Number(trackIndex);
  if (!Number.isInteger(idx) || idx < 0) return 0;
  return (Array.isArray(rows) ? rows : []).reduce(
    (count, row) => count + (depthFileSlotsFromValue(row)[idx] ? 1 : 0),
    0
  );
};

export const padDepthFileSlots = (value, width) => {
  const next = depthFileSlotsFromValue(value);
  const targetWidth = Math.max(next.length, Math.max(0, Number(width) || 0));
  while (next.length < targetWidth) next.push(null);
  return next;
};

export const clearDepthTrackSourceAt = (value, trackIndex, width = null) => {
  const idx = Number(trackIndex);
  if (!Number.isInteger(idx) || idx < 0) return depthFileSlotsFromValue(value);
  const targetWidth = width === null || width === undefined
    ? idx + 1
    : Math.max(idx + 1, Number(width) || 0);
  const next = padDepthFileSlots(value, targetWidth);
  next[idx] = null;
  return next;
};

export const removeDepthTrackColumnAt = (rows, trackIndex) => {
  const idx = Number(trackIndex);
  const sourceRows = Array.isArray(rows) ? rows : [];
  if (!Number.isInteger(idx) || idx < 0) {
    return sourceRows.map((row) => depthFileSlotsFromValue(row));
  }
  return sourceRows.map((row) => {
    const next = depthFileSlotsFromValue(row);
    if (idx < next.length) next.splice(idx, 1);
    return next;
  });
};

export const normalizeDepthTrackConfig = (entry, index, defaults = {}) => {
  const source = entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {};
  const hasHeight = Object.prototype.hasOwnProperty.call(source, 'height');
  const labelForIndex = typeof defaults.labelForIndex === 'function'
    ? defaults.labelForIndex
    : fallbackLabel;
  return {
    label: String(source.label ?? labelForIndex(index)),
    color: String(source.color || fallbackColor(index, defaults)),
    height: normalizePositiveNumber(hasHeight ? source.height : defaults.depthHeight),
    large_tick_interval: normalizePositiveNumber(
      source.large_tick_interval ?? source.tick_interval ?? defaults.largeTickInterval
    ),
    small_tick_interval: normalizePositiveNumber(
      source.small_tick_interval ?? defaults.smallTickInterval
    ),
    tick_font_size: normalizePositiveNumber(
      source.tick_font_size ?? defaults.tickFontSize
    )
  };
};

export const ensureDepthTrackConfigCount = (tracks, count, defaults = {}) => {
  const rawTracks = Array.isArray(tracks) ? tracks : [];
  const targetCount = Math.max(1, Number(count) || 1);
  const normalized = rawTracks.map((entry, index) => (
    normalizeDepthTrackConfig(entry, index, defaults)
  ));
  while (normalized.length < targetCount) {
    normalized.push(normalizeDepthTrackConfig(null, normalized.length, defaults));
  }
  return normalized;
};

export const ensureDepthTrackConfigShape = (tracks, count, defaults = {}) => {
  const targetTracks = Array.isArray(tracks) ? tracks : [];
  const targetCount = Math.max(1, Number(count) || 1);
  for (let index = 0; index < targetCount; index += 1) {
    const current = targetTracks[index];
    if (!current || typeof current !== 'object' || Array.isArray(current)) {
      targetTracks[index] = normalizeDepthTrackConfig(null, index, defaults);
      continue;
    }
    const fallback = normalizeDepthTrackConfig(null, index, defaults);
    if (current.label === undefined || current.label === null) current.label = fallback.label;
    if (current.color === undefined || current.color === null) current.color = fallback.color;
    if (current.height === undefined) current.height = fallback.height;
    if (current.large_tick_interval === undefined) current.large_tick_interval = fallback.large_tick_interval;
    if (current.small_tick_interval === undefined) current.small_tick_interval = fallback.small_tick_interval;
    if (current.tick_font_size === undefined) current.tick_font_size = fallback.tick_font_size;
  }
  return targetTracks;
};

export const reconcileDepthTracksToFiles = ({
  files,
  depthTracks,
  targetCount,
  defaults = {}
} = {}) => {
  const fileSlots = Array.isArray(files) ? files : depthFileSlotsFromValue(files);
  const rawTracks = Array.isArray(depthTracks) ? depthTracks : [];
  const usedTrackIndexes = new Set();
  return Array.from({ length: Math.max(1, Number(targetCount) || 1) }, (_, index) => {
    const file = fileSlots[index] || null;
    let sourceIndex = index;
    if (
      file &&
      rawTracks[index] &&
      !depthLabelMatchesFile(rawTracks[index]?.label, file)
    ) {
      const matchingIndex = rawTracks.findIndex((entry, candidateIndex) => (
        candidateIndex > index &&
        !usedTrackIndexes.has(candidateIndex) &&
        depthLabelMatchesFile(entry?.label, file)
      ));
      if (matchingIndex >= 0) sourceIndex = matchingIndex;
    }
    usedTrackIndexes.add(sourceIndex);
    return normalizeDepthTrackConfig(rawTracks[sourceIndex], index, defaults);
  });
};

export const removeDepthTrackAt = ({
  files,
  depthTracks,
  index,
  minCount = 1,
  defaults = {}
} = {}) => {
  const idx = Number(index);
  const originalFiles = depthFileSlotsFromValue(files);
  const originalTracks = Array.isArray(depthTracks) ? depthTracks.slice() : [];
  if (!Number.isInteger(idx) || idx < 0) {
    return {
      files: compactDepthFileSlots(originalFiles),
      depthTracks: ensureDepthTrackConfigCount(originalTracks, minCount, defaults),
      removed: false,
      removedFile: null,
      previousFiles: originalFiles
    };
  }
  const visibleCount = Math.max(minCount, originalFiles.length, 1);
  const nextFiles = originalFiles.slice();
  const nextTracks = originalTracks.slice();
  const removedFile = nextFiles[idx] || null;
  if (visibleCount <= 1) {
    nextFiles[0] = null;
    nextTracks.splice(0, nextTracks.length);
  } else {
    nextFiles.splice(idx, 1);
    if (idx < nextTracks.length) nextTracks.splice(idx, 1);
  }
  const compactedFiles = compactDepthFileSlots(nextFiles);
  const normalizedTracks = ensureDepthTrackConfigCount(
    nextTracks,
    Math.max(minCount, compactedFiles.length),
    defaults
  );
  return {
    files: compactedFiles,
    depthTracks: normalizedTracks,
    removed: true,
    removedFile,
    previousFiles: originalFiles
  };
};

export const depthSlotTrackIndex = (slot, fallbackIndex = null) => {
  const rawTrackIndex = slot?.params?.track_index;
  if (rawTrackIndex !== null && rawTrackIndex !== undefined && rawTrackIndex !== '') {
    const parsed = Number(rawTrackIndex);
    if (Number.isInteger(parsed) && parsed >= 0) return parsed;
  }
  if (slot?.depth_binding_error) return null;
  const idMatch = String(slot?.id || '').trim().match(/^depth_(\d+)$/);
  if (idMatch) return Math.max(0, Number(idMatch[1]) - 1);
  if (fallbackIndex !== null && fallbackIndex !== undefined) {
    const fallback = Number(fallbackIndex);
    if (Number.isInteger(fallback) && fallback >= 0) return fallback;
  }
  return null;
};

export const referencedDepthTrackWidth = (slots) => {
  let width = 0;
  (Array.isArray(slots) ? slots : []).forEach((slot) => {
    if (!slot || String(slot.renderer || '') !== 'depth') return;
    const trackIndex = depthSlotTrackIndex(slot);
    if (trackIndex !== null) width = Math.max(width, trackIndex + 1);
  });
  return width;
};

export const depthTrackSessionWidth = ({ rows, depthTracks, slots } = {}) => Math.max(
  depthTrackMatrixWidth(rows),
  Array.isArray(depthTracks) ? depthTracks.length : 0,
  referencedDepthTrackWidth(slots)
);

const slotHasGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.spacing) !== null ||
  normalizeOptionalText(slot?.inner_gap_px) !== null ||
  normalizeOptionalText(slot?.outer_gap_px) !== null ||
  Number(slot?.z || 0) !== 0
);

const paramsOnlyAllowedKeys = (params, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  return Object.entries(cloneParams(params)).every(([key, value]) => (
    allowed.has(String(key)) || normalizeOptionalText(value) === null
  ));
};

export const isDefaultManagedDepthSlot = (slot) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return false;
  if (String(slot.renderer || '') !== 'depth') return false;
  if (slotHasGeometry(slot)) return false;
  const id = String(slot.id || '').trim();
  if (!/^depth(?:_\d+)?$/.test(id)) return false;
  return paramsOnlyAllowedKeys(slot.params, ['track_index', 'legend_label', 'managed']);
};

const cloneSlot = (slot) => ({
  ...slot,
  params: cloneParams(slot?.params)
});

const disableInvalidManualDepthSlot = (slot, removedTrackIndex = null) => {
  const next = cloneSlot(slot);
  next.enabled = false;
  if (/^depth(?:_\d+)?$/.test(String(next.id || '').trim())) {
    next.id = `disabled_${next.id}`;
  }
  const indexText = Number.isInteger(Number(removedTrackIndex))
    ? ` ${Number(removedTrackIndex)}`
    : '';
  next.depth_binding_error = `Depth logical track index${indexText} is no longer available.`;
  delete next.params.track_index;
  delete next.params.legend_label;
  return next;
};

export const reindexDepthSlots = ({
  slots,
  removedIndex,
  activeCount,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const idx = Number(removedIndex);
  const count = Math.max(0, Number(activeCount) || 0);
  if (!Array.isArray(slots) || !Number.isInteger(idx) || idx < 0) {
    return Array.isArray(slots) ? slots.map((slot) => cloneSlot(slot)) : [];
  }
  const nextSlots = [];
  slots.forEach((slot, ordinal) => {
    if (!slot || String(slot.renderer || '') !== 'depth') {
      nextSlots.push(cloneSlot(slot));
      return;
    }
    const oldTrackIndex = depthSlotTrackIndex(slot, ordinal);
    const isManaged = managedPredicate(slot);
    if (oldTrackIndex === null || oldTrackIndex === idx) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot, oldTrackIndex ?? idx));
      return;
    }
    const newTrackIndex = oldTrackIndex > idx ? oldTrackIndex - 1 : oldTrackIndex;
    if (newTrackIndex < 0 || newTrackIndex >= count) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot, oldTrackIndex));
      return;
    }
    const next = cloneSlot(slot);
    next.params.track_index = newTrackIndex;
    nextSlots.push(next);
  });
  return nextSlots;
};

export const dropInvalidManagedDepthSlots = ({
  slots,
  activeCount,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const count = Math.max(0, Number(activeCount) || 0);
  if (!Array.isArray(slots)) return [];
  return slots
    .map((slot, ordinal) => {
      if (!slot || String(slot.renderer || '') !== 'depth') return cloneSlot(slot);
      const trackIndex = depthSlotTrackIndex(slot, ordinal);
      if (trackIndex !== null && trackIndex >= 0 && trackIndex < count) {
        const next = cloneSlot(slot);
        next.params.track_index = trackIndex;
        return next;
      }
      return managedPredicate(slot) ? null : disableInvalidManualDepthSlot(slot, trackIndex);
    })
    .filter(Boolean);
};

export const syncDepthSlotLabels = ({ slots, depthTracks, activeCount = null } = {}) => {
  if (!Array.isArray(slots)) return;
  const count = activeCount === null || activeCount === undefined
    ? (Array.isArray(depthTracks) ? depthTracks.length : 0)
    : Math.max(0, Number(activeCount) || 0);
  slots.forEach((slot, ordinal) => {
    if (!slot || String(slot.renderer || '') !== 'depth') return;
    const trackIndex = depthSlotTrackIndex(slot, ordinal);
    if (trackIndex === null || trackIndex < 0 || trackIndex >= count) {
      return;
    }
    const label = String(depthTracks?.[trackIndex]?.label ?? '').trim();
    slot.params = cloneParams(slot.params);
    slot.params.track_index = trackIndex;
    if (label) {
      slot.params.legend_label = label;
    } else {
      delete slot.params.legend_label;
    }
  });
};
